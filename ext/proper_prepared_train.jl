# Fixed scalar train lowering. The ordinary train routines remain the geometry
# authority and numerical oracle. Shared-prefix partial execution retains its
# existing executor and field/checkpoint ownership.
struct TrainPropagationSegment{P,M,T}
    propagation::P
    multipliers::M
    tail::T
end
struct TrainRuntimeOPD{I}
end

struct PreparedTrainExecution{S,G,A,P}
    segments::S
    geometry::G
    field::A
    empty::A
    phase::P
end

_train_propagation_kind(::typeof(Proper.prop_ptp)) = :ptp
_train_propagation_kind(::typeof(Proper.prop_wts)) = :wts
_train_propagation_kind(::typeof(Proper.prop_stw)) = :stw

function _prepare_train_propagations!(wf, distance, context)
    P = typeof(Proper.prepare_propagation(wf, zero(distance), context))
    propagation = P[]
    observe = function (operation, wavefront, dz, ctx)
        kind = operation === Proper.prop_stw && wavefront.reference_surface !== Proper.SPHERICAL ?
            :ptp : _train_propagation_kind(operation)
        push!(propagation, Proper.prepare_propagation(wavefront, dz, ctx; kind))
        operation(wavefront, dz, ctx)
    end
    Proper._prop_propagate!(wf, distance, Proper.PropagateOptions(false), context, observe)
    return Tuple(propagation)
end

function _prepare_train_multiplier(mask, reflection, target)
    n = size(mask, 1)
    values = Backends.host_array(mask)
    indices = Int[min(i, n - i + 2) for i in 1:n]
    quadrant = view(values, 1:fld(n, 2) + 1, 1:fld(n, 2) + 1)
    # Exact reconstruction, including signed zeros; asymmetric masks retain
    # dense storage. Reflection changes neither coordinates nor coefficients.
    if isequal(quadrant[indices, indices], values)
        compact = Backends.allocate_device_array(target, eltype(mask), size(quadrant)...)
        copyto!(compact, Array(quadrant))
        return view(compact, reflection, reflection)
    end
    return mask
end

function _prepare_train_execution(owner, target)
    AdaptiveOpticsSim.Optics._prepared_optical_train_enabled(target, owner.workspace.wavefront.field) || return nothing
    any(step -> step.op in (:propagate, :propagate_to_focus), owner.plan.steps) || return nothing
    return _prepare_train_execution(owner, target, Val(:enabled))
end

function _prepare_train_execution(owner, target, ::Val{:enabled}; cache_limit=128 * 1024 * 1024)
    cache_limit >= 0 || throw(AOG.AlgorithmGraphError("prepared train cache budget must be nonnegative"))
    plan, w = owner.plan, owner.workspace
    wf, context = w.wavefront, w.context
    # Host qualification tests may not have exercised the accelerator warmup.
    # The ordinary trajectory establishes the same final FFT alias in that case.
    if wf.field !== Proper.fft_workspace(context).scratch
        _step_optical_train_reference!(owner)
    end
    if wf.field !== Proper.fft_workspace(context).scratch
        AOG.reset_graph_node!(owner)
        return nothing
    end
    _begin_optical_train!(owner)
    n = plan.resolution
    reflection = Backends.allocate_device_array(target, Int32, n)
    copyto!(reflection, Int32[min(i, n - i + 2) for i in 1:n])
    propagation = ()
    multipliers = ()
    pending = ()
    segments = ()
    cache_bytes = sizeof(reflection)
    # Bound retained acceleration storage for other configurable trains. The
    # numerical fallback remains available when this cold resource budget is exceeded.
    for (index, step) in enumerate(plan.steps)
        if step.op === :add_opd
            segments = (segments..., TrainPropagationSegment(propagation, multipliers, pending), TrainRuntimeOPD{index}())
            propagation = ()
            multipliers = ()
            pending = ()
        elseif step.op in (:propagate, :propagate_to_focus)
            distance = step.op === :propagate ? step.distance_m : Proper.prop_get_distancetofocus(wf) + step.offset_m
            current = _prepare_train_propagations!(wf, distance, context)
            for p in current
                for phase in (p.before, p.transfer, p.after)
                    storage = Proper._propagation_phase_storage(phase)
                    storage === nothing || (cache_bytes += sizeof(storage))
                end
                p.field === wf.field || throw(AOG.AlgorithmGraphError("prepared train FFT storage differs from live field"))
                propagation = (propagation..., p)
                multipliers = (multipliers..., pending)
                pending = ()
            end
        else
            # Evaluate each fixed multiplier with the original backend and
            # geometry. Only its immutable values are cached; no optical gain
            # or adjacent multiplier is combined/reassociated during preparation.
            fill!(wf.field, one(eltype(wf.field)))
            _train_optic!(wf, step, w.masks[index], context, w.shifted_mask)
            mask = _prepare_train_multiplier(copy(wf.field), reflection, target)
            pending = (pending..., mask)
            cache_bytes += sizeof(parent(mask))
        end
        if cache_bytes > cache_limit
            AOG.reset_graph_node!(owner)
            return nothing
        end
    end
    segments = (segments..., TrainPropagationSegment(propagation, multipliers, pending))
    execution = PreparedTrainExecution(segments, _wavefront_geometry(wf), wf.field, similar(wf.field, 0, 0),
        Proper._prepare_phase(wf))
    AOG.reset_graph_node!(owner)
    return execution
end

_step_train_execution!(owner::OpticalTrainOwner, ::Nothing) = _step_optical_train_reference!(owner)
_step_train_execution!(owner::RuntimeOPDTrainOwner, ::Nothing) = _step_runtime_opd_reference!(owner)
_begin_prepared_train!(owner::OpticalTrainOwner, phase) = _begin_optical_train!(owner, phase)
_begin_prepared_train!(owner::RuntimeOPDTrainOwner, phase) = _begin_sampled_pupil!(owner.train, owner.entrance, phase)

function _step_train_execution!(owner, execution::PreparedTrainExecution)
    train = _scalar_train(owner)
    workspace = train.workspace
    workspace.wavefront.field === execution.field && Proper.fft_workspace(workspace.context).scratch === execution.field ||
        throw(AOG.AlgorithmGraphError("prepared train live field binding changed"))
    _begin_prepared_train!(owner, execution.phase)
    foreach(segment -> _step_train_segment!(owner, segment, execution.empty, execution.phase), execution.segments)
    _restore_wavefront_geometry!(train.workspace.wavefront, execution.geometry)
    _finish_optical_train!(train.workspace.wavefront, train.output, train.plan.output_quantity)
    return nothing
end
function _step_train_segment!(owner, segment::TrainPropagationSegment, empty, phase)
    train = _scalar_train(owner)
    _step_train_propagation!(train.workspace.wavefront.field, segment.propagation, segment.multipliers, segment.tail, empty)
    return nothing
end

function _step_train_propagation!(field, propagation::Tuple{}, multipliers, tail, empty)
    for mask in tail
        field .*= mask
    end
    return nothing
end

function _step_train_propagation!(field, propagation::Tuple{P,Vararg{P}}, multipliers, tail, empty) where {P}
    first(propagation).field === field || throw(AOG.AlgorithmGraphError("prepared train live field binding changed"))
    Proper._propagate_inplace!(propagation, multipliers, tail, empty)
    return nothing
end

function _step_train_segment!(owner::RuntimeOPDTrainOwner, ::TrainRuntimeOPD{I}, empty, phase) where {I}
    _step_opd_optic!(owner.train.workspace, owner.inputs[I], owner.train.plan.steps[I],
        owner.train.workspace.masks[I], owner.interpolation[I], phase)
    return nothing
end
