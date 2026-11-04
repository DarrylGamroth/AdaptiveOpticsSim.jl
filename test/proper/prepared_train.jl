# Portable numerical gate for the backend-qualified scalar train lowering.
# The ordinary executor computes its own live phases/geometry on every call.
using Test,AdaptiveOpticsSim,Proper,LinearAlgebra
using AdaptiveOpticsSim.AlgorithmGraphs

function prepared_train_fixture(array,target,T,n)
    p=4
    steps=((op="circular_aperture",radius_m=.0035),
        (op="propagate",distance_m=.04),(op="lens",focal_length_m=.08),
        (op="propagate_to_focus",offset_m=0.0),
        (op="add_opd",name="live_dm",input="dm",size=(5,6),axes=(2,1),
            sampling_m=(.00002,.000018),center_pixel=(3.2,2.7),schema="test.opd/1",outside_opd_m=0.0),
        (op="sampled_mask",name="asymmetric_phasor",amplitude="mask",phase="phase",
            sampling_m=(.00002,.000018),center_pixel=(3.2,2.7),outside_amplitude=1.0,outside_phase_rad=0.0),
        (op="lens",focal_length_m=.06),(op="propagate",distance_m=-.01),
        (op="multiply",amplitude=.93),(op="circular_obscuration",radius_m=.000001),
        (op="propagate_to_focus",offset_m=.004),(op="lens",focal_length_m=-.04),
        (op="propagate",distance_m=.03))
    config=(;resolution=n,pupil_resolution=p,diameter_m=.008,wavelength_um=1.55,steps,
        pupil_grid=(size=(p,p),axes=(2,1),sampling_m=(.0021,.0019),center_pixel=(2.7,2.4)),
        mask="mask",phase="phase")
    bindings=(mask=T[.7+.1sin(.3i+.1j)^2 for i in 1:5,j in 1:6],
        phase=T[.04i-.02j for i in 1:5,j in 1:6])
    inputs=(opd=array(zeros(T,p,p)),amplitude=array(ones(T,p,p)),dm=array(zeros(T,5,6)))
    function prepare()
        node=proper_optical_train_node_factory(;T)(:optics,config,(;),bindings)
        definition=algorithm_graph((node,);inputs=(graph_input(:opd,:optics=>:pupil_opd,inputs.opd),
            graph_input(:amplitude,:optics=>:pupil_amplitude,inputs.amplitude),
            graph_input(:dm,:optics=>:dm,inputs.dm)),outputs=(graph_output(:intensity,:optics=>:output),))
        prepare_algorithm_graph(definition;target,execution=StreamGraphExecution())
    end
    return (;inputs,prepare)
end

function force_prepared_train(owner,target)
    ext=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimProperExt)
    train=ext._scalar_train(owner);w=train.workspace
    execution=ext._prepare_train_execution(train,target,Val(:enabled))
    execution===nothing && error("fixture did not establish raw FFT storage")
    workspace=ext.OpticalTrainWorkspace(w.wavefront,w.context,w.padded_opd,w.padded_amplitude,
        w.masks,w.shifted_mask,execution)
    replacement=ext.OpticalTrainOwner(train.plan,workspace,train.pupil_opd,train.pupil_amplitude,
        train.output,train.sampling_m)
    return ext.RuntimeOPDTrainOwner(replacement,owner.interpolation,owner.inputs,owner.entrance)
end

function check_prepared_train(array,target;captured=false)
    ext=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimProperExt)
    backend=AdaptiveOpticsSim.Backends
    @testset "Prepared scalar train $array" begin
        for T in (Float32,Float64),n in (8,16)
            fixture=prepared_train_fixture(array,target,T,n)
            reference_graph=fixture.prepare();candidate_graph=fixture.prepare()
            reference=AlgorithmGraphs.prepared_graph_node(reference_graph,:optics)
            owner=backend._with_prepared_device_execution_context(candidate_graph.context) do
                force_prepared_train(AlgorithmGraphs.prepared_graph_node(candidate_graph,:optics),target)
            end
            field=owner.train.workspace.wavefront.field
            operation()=AlgorithmGraphs.step_graph_node!(owner)
            recorded=backend._with_prepared_device_execution_context(candidate_graph.context) do
                operation();backend._synchronize_prepared_device_execution_context!(candidate_graph.context)
                captured ? backend._capture_prepared_device_graph(operation,candidate_graph.context) : nothing
            end
            first_images=nothing
            for state in (1,2,1)
                copyto!(fixture.inputs.opd,array(T[state*2e-8sin(.2i+.3j) for i in 1:4,j in 1:4]))
                copyto!(fixture.inputs.amplitude,array(T[.6+.3cos(state*(.1i+.2j))^2 for i in 1:4,j in 1:4]))
                copyto!(fixture.inputs.dm,array(T[state*3e-8cos(.2i-.1j) for i in 1:5,j in 1:6]))
                backend._synchronize_prepared_device_execution_context!(candidate_graph.context)
                backend._with_prepared_device_execution_context(reference_graph.context) do
                    fill!(reference.train.workspace.wavefront.field,Complex{T}(NaN))
                    ext._step_runtime_opd_reference!(reference)
                    backend._synchronize_prepared_device_execution_context!(reference_graph.context)
                end
                backend._with_prepared_device_execution_context(candidate_graph.context) do
                    fill!(field,Complex{T}(NaN))
                    recorded===nothing ? operation() : backend._launch_prepared_device_graph!(recorded,candidate_graph.context)
                    backend._synchronize_prepared_device_execution_context!(candidate_graph.context)
                end
                expected=Array(reference.train.output);actual=Array(owner.train.output)
                tolerance=T===Float32 ? 1e-5 : 1e-11
                @test norm(actual-expected)/norm(expected)<=tolerance
                @test maximum(abs,actual-expected)/maximum(abs,expected)<=tolerance
                expected_field=Array(reference.train.workspace.wavefront.field)
                @test norm(Array(field)-expected_field)/norm(expected_field)<=tolerance
                @test ext._wavefront_geometry(owner.train.workspace.wavefront)==ext._wavefront_geometry(reference.train.workspace.wavefront)
                @test owner.train.workspace.wavefront.field===field
                first_images===nothing && (first_images=copy(actual))
                state==1 && @test actual==first_images
            end
            # A rejected cache budget preserves the ordinary fallback contract.
            disabled=backend._with_prepared_device_execution_context(reference_graph.context) do
                ext._prepare_train_execution(reference.train,target,Val(:enabled);cache_limit=0)
            end
            @test disabled===nothing
            backend._with_prepared_device_execution_context(reference_graph.context) do
                ext._step_runtime_opd_reference!(reference)
                backend._synchronize_prepared_device_execution_context!(reference_graph.context)
            end
            @test Array(reference.train.output)≈first_images rtol=(T===Float32 ? 1e-5 : 1e-11)
            # Field bindings cannot drift silently into an existing capture.
            owner.train.workspace.wavefront.field=copy(field)
            @test_throws AlgorithmGraphError AlgorithmGraphs.step_graph_node!(owner)
            owner.train.workspace.wavefront.field=field
            AlgorithmGraphs.reset_graph_node!(owner)
            @test all(iszero,Array(owner.train.output))
            @test owner.train.workspace.wavefront.field===field
            backend._with_prepared_device_execution_context(candidate_graph.context) do
                recorded===nothing ? operation() : backend._launch_prepared_device_graph!(recorded,candidate_graph.context)
                backend._synchronize_prepared_device_execution_context!(candidate_graph.context)
            end
            @test Array(owner.train.output)==first_images
        end
    end
end
