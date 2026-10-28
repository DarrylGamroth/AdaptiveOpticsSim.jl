# Run after test_zemax.jl in the zemax selector; reuse its authored source builder.
@testset "authored static optical variants preserve base graph bytes" begin
    mktempdir() do root
        write_zemax_fixture(root)
        Z.prepare(root;yamlfmt="true")
        base=read(joinpath(root,"generated/test.yaml"))
        path=joinpath(root,"simulation/definition.yaml")
        definition=YAML.load_file(path)
        circles=[Dict("radius_m"=>0.0005,"x_m"=>0.0,"y_m"=>0.0)]
        definition["preparation"]["variants"]=Dict("test_unfringed"=>Dict("branch"=>"test",
            "settings"=>Dict("masks"=>Dict("lyot"=>Dict("circles"=>circles)))))
        YAML.write_file(path,definition)
        authored=authored_fixture_bytes(root)
        @test Z.prepare(root;yamlfmt="true")===nothing
        @test read(joinpath(root,"generated/test.yaml"))==base
        @test authored_fixture_bytes(root)==authored
        runtime=YAML.load_file(joinpath(root,"generated/simulation.yaml"))
        @test runtime["graphs"]["test_unfringed"]=="generated/test_unfringed.yaml"
        @test runtime["optical_model"]["variants"]==Dict("test_unfringed"=>Dict("branch"=>"test"))
        @test !haskey(runtime,"preparation")
        variant=YAML.load_file(joinpath(root,"generated/test_unfringed.yaml"))
        @test variant["name"]=="synthetic_test_unfringed"
        mask=only(step for step in only(variant["nodes"])["config"]["steps"] if step["op"]=="circular_mask")
        @test mask["circles"]==circles && mask["mode"]=="transmit"
        receipt=Z.validate_provenance(root)
        @test Set(entry["path"] for entry in receipt["graphs"])==Set(["generated/test.yaml","generated/test_unfringed.yaml"])
        generated=Dict(relative=>read(joinpath(root,relative)) for relative in
            ("generated/test.yaml","generated/test_unfringed.yaml","generated/simulation.yaml","generated/optics/provenance.yaml"))
        replanned=Z.execution_graphs(root)
        @test replanned["test_unfringed"]==variant
        @test bounded_zemax_storage(replanned)
        @test Z.prepare(root;yamlfmt="true")===nothing
        @test all(read(joinpath(root,relative))==bytes for (relative,bytes) in generated)
        write(joinpath(root,"generated/test_unfringed.yaml"),"tampered")
        @test_throws ArgumentError Z.validate_provenance(root)
        @test Z.prepare(root;yamlfmt="true")===nothing
        @test Z.validate_provenance(root)["schema_version"]==1
        # Updating receipt hashes cannot authorize a forged physical branch.
        runtime_path=joinpath(root,"generated/simulation.yaml")
        forged=YAML.load_file(runtime_path)
        forged["optical_model"]["variants"]["test_unfringed"]["branch"]="other"
        YAML.write_file(runtime_path,forged)
        forged_receipt=YAML.load_file(joinpath(root,"generated/optics/provenance.yaml"))
        forged_receipt["execution"]["sha256"]=bytes2hex(sha256(read(runtime_path)))
        @test_throws ArgumentError Z.validate_provenance(root,forged_receipt)
        @test Z.prepare(root;yamlfmt="true")===nothing
        @test Z.validate_provenance(root)["schema_version"]==1
    end
end

@testset "static variants retain physical spectral and filter routing" begin
    mktempdir() do root
        write_zemax_fixture(root)
        path=joinpath(root,"simulation/definition.yaml")
        definition=YAML.load_file(path)
        definition["defaults"]["spectra"]=Dict("test"=>Dict("wavelength_um"=>1.55,
            "fractional_bandwidth"=>0.02,"samples"=>3,"filter"=>"attenuator"))
        definition["preparation"]["variants"]=Dict("test_unfringed"=>Dict("branch"=>"test",
            "settings"=>Dict("masks"=>Dict("lyot"=>Dict("circles"=>[
                Dict("radius_m"=>0.0005,"x_m"=>0.0,"y_m"=>0.0)])))))
        YAML.write_file(path,definition)
        change_physical_fixture(root, instrument -> begin
            instrument["optical_filters"]=Dict("attenuator"=>Dict(
                "wavelength_m"=>[1.0e-6,2.0e-6],"transmission"=>[0.25,0.25]))
            instrument["optical_paths"]["test"]["filters"]=["attenuator"]
        end)
        Z.prepare(root;yamlfmt="true")
        before=authored_fixture_bytes(root)
        snapshots=Dict(p=>read(joinpath(root,p)) for p in
            ("generated/test.yaml","generated/test_unfringed.yaml","generated/simulation.yaml","generated/optics/provenance.yaml"))
        settings=Dict("resolution"=>32,"spectra"=>Dict("test"=>Dict("wavelength_um"=>1.7)))
        graphs=Z.execution_graphs(root;settings)
        @test Set(keys(graphs))==Set(["test","test_unfringed"])
        base=only(graphs["test"]["nodes"])["config"]
        variant=only(graphs["test_unfringed"]["nodes"])["config"]
        @test base["resolution"]==variant["resolution"]==32
        @test [c["wavelength_um"] for c in base["channels"]]==[c["wavelength_um"] for c in variant["channels"]]
        @test all(c["transmission"]==0.25 for c in variant["channels"])
        for channel in variant["channels"]
            mask=only(s for s in channel["steps"] if s["op"]=="circular_mask")
            @test only(mask["circles"])["radius_m"]==0.0005 && mask["mode"]=="transmit"
        end
        @test authored_fixture_bytes(root)==before
        @test all(read(joinpath(root,p))==bytes for (p,bytes) in snapshots)
        @test bounded_zemax_storage(graphs)
    end
end

@testset "invalid static variant authoring cannot replace generated files" begin
    variants=[
        [],
        Dict("a$i"=>Dict("branch"=>"test","settings"=>Dict()) for i in 1:17),
        Dict("../outside"=>Dict("branch"=>"test","settings"=>Dict())),
        Dict("Upper"=>Dict("branch"=>"test","settings"=>Dict())),
        Dict(repeat("a",65)=>Dict("branch"=>"test","settings"=>Dict())),
        Dict("test"=>Dict("branch"=>"test","settings"=>Dict())),
        Dict("simulation"=>Dict("branch"=>"test","settings"=>Dict())),
        Dict("variant"=>Dict("branch"=>"missing","settings"=>Dict())),
        Dict("variant"=>Dict("branch"=>"test")),
        Dict("variant"=>Dict("branch"=>"test","settings"=>Dict(),"unknown"=>1)),
        Dict("variant"=>Dict("branch"=>"test","settings"=>1)),
        Dict("variant"=>Dict("branch"=>"test","settings"=>Dict("unknown"=>1))),
        Dict("variant"=>Dict("branch"=>"test","settings"=>Dict("masks"=>Dict("lyot"=>Dict("unknown"=>1)))))
    ]
    mktempdir() do root
        write_zemax_fixture(root)
        Z.prepare(root;yamlfmt="true")
        path=joinpath(root,"simulation/definition.yaml")
        definition=YAML.load_file(path)
        original=Dict(p=>read(joinpath(root,p)) for p in
            ("generated/test.yaml","generated/simulation.yaml","generated/optics/provenance.yaml"))
        for value in variants
            definition["preparation"]["variants"]=value
            YAML.write_file(path,definition)
            @test_throws ArgumentError Z.prepare(root;yamlfmt="true")
            @test all(read(joinpath(root,p))==bytes for (p,bytes) in original)
        end
        definition["preparation"]["variants"]=Dict("a$i"=>Dict("branch"=>"test","settings"=>Dict()) for i in 1:16)
        YAML.write_file(path,definition)
        @test Z.prepare(root;yamlfmt="true")===nothing
        @test length(Z.validate_provenance(root)["graphs"])==17
        @test read(joinpath(root,"generated/test.yaml"))==original["generated/test.yaml"]
    end
end
