# Offline regression against hash-bound outputs from the original cold compiler.
# Agreement checks shared ancestry; it does not establish physical accuracy.
using AdaptiveOpticsSim,Proper,FITSIO,YAML,SHA,TOML,LinearAlgebra,Test
using AdaptiveOpticsSim.AlgorithmGraphs
length(ARGS) in (2,3) || error("Usage: compare_compiled_optical_paths.jl MODEL.yaml FROZEN_EXPORT [Float32|Float64]")
model_path,legacy_root=abspath.(ARGS[1:2])
T=length(ARGS)==3 ? Dict("Float32"=>Float32,"Float64"=>Float64)[ARGS[3]] : Float64
hash(path)=bytes2hex(sha256(read(path)))
metadata,reader=AlgorithmGraphs._read_compiled_optical_model(model_path,hash(model_path))
legacy_manifest=joinpath(legacy_root,"manifest.toml")
metadata["provenance"]["imported_manifest_sha256"]==hash(legacy_manifest) || error("reference manifest identity differs")
legacy=TOML.parsefile(legacy_manifest)
function artifact(reference,expected)
    !isabspath(reference) && normpath(reference)==reference && !(".." in splitpath(reference)) || error("invalid reference path")
    path=joinpath(legacy_root,reference)
    hash(path)==expected || error("reference artifact hash differs")
    return FITS(path,"r") do f;permutedims(read(f[1]),(2,1));end
end
live=metadata["live"]
dm=zeros(T,live["dm_size"],live["dm_size"])
opd=zeros(T,live["entrance_size"],live["entrance_size"])
a=zeros(T,size(opd))
owner=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(metadata,reader;
    T,dm_opd=dm,entrance_opd=opd,entrance_attenuation=a)
tolerance=T===Float32 ? 5e-4 : 2e-8
@testset "Compiled optical model frozen reference regression" begin
    for scenario in legacy["references"]
        dm.=artifact(scenario["dm_opd"],scenario["dm_opd_sha256"])
        if haskey(scenario,"entrance_opd")
            opd.=artifact(scenario["entrance_opd"],scenario["entrance_opd_sha256"])
            a.=artifact(scenario["entrance_attenuation"],scenario["entrance_attenuation_sha256"])
        else
            fill!(opd,0);fill!(a,0)
        end
        products=AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
        for name in (:fringed,:unfringed,:reference)
            key=String(name)
            expected=artifact(scenario[key],scenario[key*"_sha256"])
            actual=getproperty(products,name)
            relative=norm(actual-expected)/norm(expected)
            println(scenario["label"]," ",name," relative_error=",relative)
            @test all(isfinite,actual)
            @test relative<tolerance
        end
    end
end
