# julia --project=BENCH_ENV fixture.jl OUTPUT_DIRECTORY
using AdaptiveOpticsSim, REVOLTCopperSim
using AdaptiveOpticsSim.Optics, AdaptiveOpticsSim.WavefrontSensors
using Random, LinearAlgebra, TOML, SHA
const FixtureWFS = AdaptiveOpticsSim.WavefrontSensors
function save_array(directory, name, array)
    path=joinpath(directory,name*".f32")
    values=Float32.(array)
    open(io->write(io,vec(values)),path,"w")
    return Dict("file"=>basename(path),"shape"=>collect(size(values)),
        "dtype"=>"float32","order"=>"F","sha256"=>bytes2hex(sha256(read(path))))
end
function create_fixture(directory)
    isdir(directory) && !isempty(readdir(directory)) && error("fixture output must be empty")
    mkpath(directory)
    n=240
    tel=Telescope(resolution=n,diameter=1.22,central_obstruction=0.,T=Float32)
    src=Source(band=:V,magnitude=5.752574989159953,T=Float32)
    pupil=PupilFunction(tel;T=Float32)
    amplitude=Array(pupil.amplitude)
    coords=REVOLTCopperSim.actuator_coordinates()
    topology=AdaptiveOpticsSim.Optics.SampledActuatorTopology(coords;T=Float32)
    dm=DeformableMirror(tel;topology,influence_width=0.08626550214129701,T=Float32)
    basis=reshape(AdaptiveOpticsSim.Optics.materialize_influence_matrix(dm),n,n,277)
    names=["flat","piston_100nm","tilt_x_20nm","tilt_y_20nm","focus_20nm",
        "astigmatism_20nm","coma_20nm","tilt_x_100nm","focus_100nm",
        "fourier_4_20nm","fourier_8_20nm","fourier_12_20nm","mixed_100nm",
        "bandlimited_100nm_a","bandlimited_100nm_b","bandlimited_300nm_a","bandlimited_300nm_b"]
    probes=zeros(Float32,n,n,length(names))
    axis=(2 .* collect(1:n) .- n .- 1)./n
    for j in 1:n,i in 1:n
        x,y=axis[i],axis[j]
        probes[i,j,2]=100e-9
        probes[i,j,3]=20e-9*x;probes[i,j,4]=20e-9*y
        probes[i,j,5]=20e-9*(2*(x*x+y*y)-1)
        probes[i,j,6]=20e-9*(x*x-y*y)
        probes[i,j,7]=20e-9*(3*(x*x+y*y)-2)*x
        probes[i,j,8]=100e-9*x;probes[i,j,9]=100e-9*(2*(x*x+y*y)-1)
        for (k,fx,fy) in ((10,4,3),(11,8,5),(12,12,7))
            probes[i,j,k]=20e-9*sinpi(fx*x)*cospi(fy*y)
        end
        probes[i,j,13]=100e-9*(.5sinpi(x+2y)+.3cospi(3x-y)+.2sinpi(7x+5y))
    end
    rng=Xoshiro(823)
    support=amplitude .> 0
    for k in 14:17
        field=zeros(Float64,n,n)
        for _ in 1:48
            fx,fy=rand(rng,1:12),rand(rng,1:12)
            phase=2pi*rand(rng)
            weight=randn(rng)/(fx^2+fy^2)^(11/12)
            for j in 1:n,i in 1:n
                field[i,j]+=weight*sinpi(fx*axis[i]+fy*axis[j]+phase/pi)
            end
        end
        field .-= sum(field[support])/count(support)
        rms=sqrt(sum(abs2,field[support])/count(support))
        probes[:,:,k].=field.*((k<16 ? 100e-9 : 300e-9)/rms)
    end
    contract=Dict("resolution"=>240,"modulation_points"=>96,
        "modulation_radius_lambda_over_d"=>5.0,"diameter_m"=>1.22,"wavelength_m"=>550e-9,
        "pupil_samples"=>30,"separation_pixels"=>2,"edge_pixels"=>1,
        "photon_irradiance"=>Float64(FixtureWFS.photon_irradiance(src)),
        "exposure_seconds"=>.002,"fft_side"=>512,"output_side"=>64,
        "precision"=>"float32","calibration_poke_m"=>20e-9,
        "calibration_columns"=>277,"input_axes"=>"axis0=x, axis1=y; identical array ordering across engines",
        "scope"=>"Actual Copper optical prescription; common supplied pupil and provisional Gaussian277 OPD modes. No detector, atmospheric propagation, DM runtime, RTC, or fitted photometric/registration correction.")
    result=Dict("contract"=>contract,"probe_names"=>names,
        "amplitude"=>save_array(directory,"amplitude",amplitude),
        "probes"=>save_array(directory,"probes",probes),
        "calibration_basis"=>save_array(directory,"calibration_basis",basis),
        "actuator_coordinates"=>save_array(directory,"actuator_coordinates",coords),
        "generator_sha256"=>bytes2hex(sha256(read(@__FILE__))),"seed"=>823,
        "aos_revision"=>readchomp(`git -C $(dirname(dirname(pathof(AdaptiveOpticsSim)))) rev-parse HEAD`),
        "copper_revision"=>readchomp(`git -C $(dirname(dirname(pathof(REVOLTCopperSim)))) rev-parse HEAD`),
        "calibration_basis_scope"=>"AOS sampled provisional Copper Gaussian operator, unit surface OPD coefficients. Shared as input; external DM implementations are not tested.",
        "bandlimited_probe_scope"=>"Deterministic finite Fourier probes, not a validated atmospheric realization.")
    open(io->TOML.print(io,result;sorted=true),joinpath(directory,"fixture.toml"),"w")
    return result
end
if abspath(PROGRAM_FILE)==@__FILE__
    create_fixture(only(ARGS))
end
