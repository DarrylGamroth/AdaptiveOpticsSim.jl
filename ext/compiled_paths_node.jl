# Generic artifact adapter for the lowered multi-path Fresnel executor. Its live
# input grids are independent of propagation and fitted-screen sample counts.
using FixedSizeArrays: FixedSizeVector
struct CompiledPathsNode{T} end
struct CompiledPathsConfiguration{D,P,R,T,C}
    path::String
    sha256::String
    dm_grid::D
    pupil_grid::P
    reference::R
    schemas::NTuple{4,String}
    mode::Symbol
    amplitude_mode::Symbol
    summation::Symbol
    output_size::Int
    sampling_m::T
    output_sampling_calibration::C
end
function AOG._compiled_optical_paths_node(::Type{T},name,config,props,bindings) where {T<:Union{Float32,Float64}}
    AOG.graph_node_bindings(AOG.CompiledOpticalPathsNodeFactory{T}(),config)
    AOG._require_named_fields(props,(),(),"compiled optical paths props")
    path=AOG._file_string(config.model_path,"compiled optical model path")
    hash=AOG._file_string(config.model_sha256,"compiled optical model SHA-256")
    metadata,_=AOG._read_compiled_optical_model(path,hash)
    dm_grid=_entrance_pupil_grid(config.dm_grid,T,AOG._file_integer(first(config.dm_grid.size),"DM resolution"))
    pupil_grid=_entrance_pupil_grid(config.pupil_grid,T,AOG._file_integer(first(config.pupil_grid.size),"pupil resolution"))
    mode=Symbol(metadata["mode"])
    mode in (:sky,:source) || throw(AOG.AlgorithmGraphError("compiled model entrance mode must be sky or source"))
    amplitude_mode=hasproperty(config,:entrance_amplitude_mode) ? Symbol(config.entrance_amplitude_mode) : :reference_ratio
    reference=amplitude_mode===:reference_ratio ? _mask_snapshot(AOG._require_binding(bindings,
        config.entrance_reference_amplitude,"entrance reference amplitude"),T,:amplitude) : nothing
    reference===nothing || size(reference)==pupil_grid.size || throw(DimensionMismatch("entrance reference amplitude grid differs"))
    summation=hasproperty(config,:summation) ? Symbol(config.summation) : :path_pairs
    summation in (:path_pairs,:coherent_modes) || throw(AOG.AlgorithmGraphError("unsupported compiled optical summation"))
    schemas=map((:dm_opd_schema,:pupil_opd_schema,:pupil_amplitude_schema,:output_schema),
        ("org.adaptiveopticssim.dm-opd/1",DEFAULT_OPD_SCHEMA,DEFAULT_AMPLITUDE_SCHEMA,DEFAULT_INTENSITY_SCHEMA)) do key,default
        hasproperty(config,key) ? AOG._file_string(getproperty(config,key),"compiled paths $key") : default
    end
    calibration=AOG._compiled_output_sampling_calibration(config,T)
    sampling=calibration===nothing ? T(metadata["dx_output"]) : calibration.sampling_m
    c=CompiledPathsConfiguration(path,hash,dm_grid,pupil_grid,reference,schemas,mode,amplitude_mode,
        summation,Int(metadata["image_size"]),sampling,calibration)
    return AOG.algorithm_node(name,CompiledPathsNode{T},c;props=NamedTuple())
end
function AOG.graph_node_ports(::Type{CompiledPathsNode{T}},c::CompiledPathsConfiguration) where {T}
    outputs=map((:fringed,:unfringed,:reference)) do role
        AOG.graph_port_contract(role,:output,:data,T,(c.output_size,c.output_size),c.schemas[4],:column_major)
    end
    return (AOG.graph_port_contract(:dm_opd,:input,:data,T,c.dm_grid.size,c.schemas[1],:column_major),
        AOG.graph_port_contract(:pupil_opd,:input,:data,T,c.pupil_grid.size,c.schemas[2],:column_major),
        AOG.graph_port_contract(:pupil_amplitude,:input,:data,T,c.pupil_grid.size,c.schemas[3],:column_major),outputs...)
end
struct CompiledPathsInputValidation{D,H}
    device::D
    host::H
end
struct CompiledPathsInputWorkspace{A,V}
    amplitude::A
    reference::A
    validation::V
end
struct CompiledPathsOwner{P,I,O,M,W,R,T,C}
    propagation::P
    inputs::I
    outputs::O
    mapping::M
    workspace::W
    reference::R
    sampling_m::T
    output_sampling_calibration::C
end
function _compiled_paths_ratio!(attenuation,amplitude,reference)
    for i in eachindex(attenuation,amplitude,reference)
        a,r=amplitude[i],reference[i]
        if iszero(r)
            iszero(a) || throw(ArgumentError("nonzero entrance amplitude outside its declared reference support"))
            attenuation[i]=zero(a)
        elseif iszero(a)
            throw(ArgumentError("new opaque entrance samples inside reference support are unsupported; regenerate the optical model for changed pupil support"))
        else
            logratio=log(a)-log(r)
            isfinite(exp(logratio)) || throw(ArgumentError("entrance amplitude ratio is not representable"))
            attenuation[i]=-logratio
        end
    end
    return nothing
end
function AOG.prepare_graph_node(::Type{CompiledPathsNode{T}},c::CompiledPathsConfiguration,
        ::NamedTuple{()},inputs::NamedTuple,outputs::NamedTuple,::NamedTuple{()},target) where {T}
    _require_compiled_input_target(target)
    for (array,shape,role) in ((inputs.dm_opd,c.dm_grid.size,"DM OPD"),
            (inputs.pupil_opd,c.pupil_grid.size,"entrance OPD"),
            (inputs.pupil_amplitude,c.pupil_grid.size,"entrance amplitude"),
            ((getproperty(outputs,name),(c.output_size,c.output_size),String(name))
                for name in (:fringed,:unfringed,:reference))...)
        Base.has_offset_axes(array) && throw(AOG.AlgorithmGraphError("compiled optical $role requires one-based storage"))
        _validate_compiled_paths_buffer(target,array,T,shape,"compiled optical $role")
    end
    metadata,read_array=AOG._read_compiled_optical_model(c.path,c.sha256)
    live=metadata["live"]
    ndm=Int(live["dm_size"])
    nentrance=Int(live["entrance_size"])
    dm=_compiled_paths_array(target,T,ndm,ndm)
    entrance=_compiled_paths_array(target,T,nentrance,nentrance)
    attenuation=_compiled_paths_array(target,T,nentrance,nentrance)
    propagation=_prepare_compiled_paths(metadata,read_array,nothing;T,target,dm_opd=dm,
        entrance_opd=entrance,entrance_attenuation=attenuation,summation=c.summation)
    dm_mapping=_prepare_grid_interpolation(
        _bilinear_interpolation_plan(c.dm_grid,ndm,T(live["dm_pitch"]),zero(T)),target)
    pupil_mapping=c.mode===:sky ? _prepare_grid_interpolation(
        _bilinear_interpolation_plan(c.pupil_grid,nentrance,T(live["entrance_pitch"]),zero(T)),target) : nothing
    reference=_compiled_paths_input_reference(c.reference,target,T)
    native_reference=_compiled_paths_array(target,T,nentrance,nentrance)
    native_amplitude=_compiled_paths_array(target,T,nentrance,nentrance)
    if c.mode===:sky && reference!==nothing
        _interpolate_grid!(native_reference,reference,pupil_mapping)
    else
        fill!(native_reference,1)
    end
    owner=CompiledPathsOwner(propagation,inputs,outputs,(dm=dm_mapping,pupil=pupil_mapping),
        CompiledPathsInputWorkspace(native_amplitude,native_reference,
            _compiled_paths_input_validation(target)),reference,c.sampling_m,c.output_sampling_calibration)
    return owner
end
_compiled_paths_input_reference(::Nothing,target,T)=nothing
_compiled_paths_input_reference(reference,target,T)=
    _compiled_paths_resident(target,Matrix(reference),T)

_compiled_paths_input_validation(::Backends.HostComputeDevice)=nothing
function _compiled_paths_input_validation(target)
    device=_compiled_paths_array(target,UInt32,6)
    # Vector backing supports the backend's ordinary bulk device-to-host copy.
    host=FixedSizeVector{UInt32,Vector{UInt32}}(undef,6)
    fill!(host,UInt32(0))
    return CompiledPathsInputValidation(device,host)
end

function _enqueue_compiled_paths_node!(owner)
    p=owner.propagation
    _interpolate_grid!(p.inputs.dm_opd,owner.inputs.dm_opd,owner.mapping.dm)
    _compiled_paths_entrance!(owner,owner.mapping.pupil)
    products=_enqueue_compiled_paths!(p)
    for name in (:fringed,:unfringed,:reference)
        Backends.copyto_backend_async!(getproperty(owner.outputs,name),getproperty(products,name))
    end
    return nothing
end
function _step_compiled_paths_node!(::Backends.ScalarCPUStyle,owner)
    _compiled_paths_validate_inputs!(owner,owner.reference)
    return _enqueue_compiled_paths_node!(owner)
end
function _step_compiled_paths_node!(style::Backends.AcceleratorStyle,owner)
    _enqueue_compiled_paths_inputs!(style,owner)
    _enqueue_compiled_paths_node!(owner)
    _check_compiled_paths_inputs!(owner.workspace.validation)
    return nothing
end
AOG.step_graph_node!(owner::CompiledPathsOwner)=
    _step_compiled_paths_node!(Backends.execution_style(owner.inputs.dm_opd),owner)
function AOG.enqueue_captured_graph_node!(owner::CompiledPathsOwner)
    _enqueue_compiled_paths_inputs!(Backends.execution_style(owner.inputs.dm_opd),owner)
    return _enqueue_compiled_paths_node!(owner)
end
_enqueue_captured_scc_product!(owner::CompiledPathsOwner, product::Val) =
    _enqueue_captured_scc_product!(owner, (product, nothing))

function _enqueue_captured_scc_product!(owner::CompiledPathsOwner, selection::Tuple{Val{P},M}) where {P,M}
    product, merge = selection
    _enqueue_compiled_paths_inputs!(Backends.execution_style(owner.inputs.dm_opd), owner)
    propagation = owner.propagation
    _interpolate_grid!(propagation.inputs.dm_opd, owner.inputs.dm_opd, owner.mapping.dm)
    _compiled_paths_entrance!(owner, owner.mapping.pupil)
    intensity = merge === nothing ? _compiled_paths_selected!(propagation, product) :
        _compiled_paths_merged!(propagation, product, merge)
    Backends.copyto_backend_async!(getproperty(owner.outputs, P), intensity)
    return nothing
end

AOG._complete_captured_graph_node!(owner::CompiledPathsOwner)=
    _check_compiled_paths_inputs!(owner.workspace.validation)
_compiled_paths_entrance!(owner,::Nothing)=nothing
function _compiled_paths_entrance!(owner,mapping::BilinearInterpolationPlan)
    inputs,p=owner.inputs,owner.propagation
    _interpolate_grid!(p.inputs.entrance_opd,inputs.pupil_opd,mapping)
    _compiled_paths_amplitude!(owner,owner.reference,mapping)
    return nothing
end
function _compiled_paths_amplitude!(owner,::Nothing,mapping)
    a=owner.inputs.pupil_amplitude
    _interpolate_grid!(owner.propagation.inputs.entrance_attenuation,a,mapping)
    return nothing
end
function _compiled_paths_amplitude!(owner,reference::AbstractMatrix,mapping)
    a=owner.inputs.pupil_amplitude
    _interpolate_grid!(owner.workspace.amplitude,a,mapping)
    _compiled_paths_ratio!(Backends.execution_style(owner.workspace.amplitude),
        owner.propagation.inputs.entrance_attenuation,owner.workspace.amplitude,
        owner.workspace.reference,owner.workspace.validation)
    return nothing
end
function AOG.reset_graph_node!(owner::CompiledPathsOwner)
    for name in (:fringed,:unfringed,:reference)
        fill!(getproperty(owner.outputs,name),0)
    end
    _reset_compiled_paths_validation!(owner.workspace.validation)
    return nothing
end
AOG.optical_sampling_m(owner::CompiledPathsOwner)=owner.sampling_m
_compiled_paths_capture_capability(target)=AOG.GraphNodeCaptureUnsupported()
_compiled_paths_capture_capability(::Backends.AcceleratorComputeDevice{<:Union{Backends.CUDABackend,Backends.AMDGPUBackend}})=
    AOG.GraphNodeCaptureSafe()
AOG.graph_node_capture_capability(owner::CompiledPathsOwner)=
    _compiled_paths_capture_capability(Backends.compute_device(owner.inputs.dm_opd))

function _compiled_paths_validate_inputs!(owner,reference)
    all(isfinite,owner.inputs.dm_opd) || throw(ArgumentError("DM OPD must be finite"))
    all(isfinite,owner.inputs.pupil_opd) || throw(ArgumentError("entrance OPD must be finite"))
    _compiled_paths_validate_amplitude!(owner.inputs.pupil_amplitude,reference)
    return nothing
end
function _compiled_paths_validate_amplitude!(amplitude,::Nothing)
    all(isfinite,amplitude) || throw(ArgumentError("entrance log attenuation must be finite"))
    return nothing
end
function _compiled_paths_validate_amplitude!(amplitude,reference::AbstractMatrix)
    for i in eachindex(amplitude,reference)
        a,r=amplitude[i],reference[i]
        isfinite(a) && a>=0 || throw(ArgumentError("entrance field amplitude must be finite and nonnegative"))
        if iszero(r)
            iszero(a) || throw(ArgumentError("nonzero entrance amplitude where its reference is opaque"))
        elseif iszero(a)
            throw(ArgumentError("new opaque entrance samples inside reference support are unsupported; regenerate the optical model for changed pupil support"))
        end
    end
    return nothing
end

_require_compiled_input_target(::Backends.HostComputeDevice)=nothing
_require_compiled_input_target(::Backends.AcceleratorComputeDevice{<:Union{Backends.CUDABackend,Backends.AMDGPUBackend}})=nothing
_require_compiled_input_target(target)=throw(AOG.AlgorithmGraphError("compiled optical input-grid adapter requires CPU, CUDA or AMDGPU execution"))

KernelAbstractions.@kernel function _compiled_paths_validate_finite_kernel!(failure,input,slot)
    index=@index(Global,Linear)
    if index <= length(input) && !isfinite(input[index])
        KernelAbstractions.@atomic failure[slot]=UInt32(1)
    end
end
@inline _compiled_paths_amplitude_failure(amplitude,::Nothing,index)=
    isfinite(amplitude) ? 0 : 3
@inline function _compiled_paths_amplitude_failure(amplitude,reference,index)
    isfinite(amplitude) && amplitude >= 0 || return 3
    r=reference[index]
    iszero(r) && !iszero(amplitude) && return 4
    !iszero(r) && iszero(amplitude) && return 5
    return 0
end
KernelAbstractions.@kernel function _compiled_paths_validate_amplitude_kernel!(failure,amplitude,reference)
    index=@index(Global,Linear)
    if index <= length(amplitude)
        slot=_compiled_paths_amplitude_failure(amplitude[index],reference,index)
        if !iszero(slot)
            KernelAbstractions.@atomic failure[slot]=UInt32(1)
        end
    end
end
function _enqueue_compiled_paths_inputs!(style,owner)
    failure=owner.workspace.validation.device
    fill!(failure,UInt32(0))
    for (input,slot) in ((owner.inputs.dm_opd,1),(owner.inputs.pupil_opd,2))
        Backends.launch_kernel_async!(style,_compiled_paths_validate_finite_kernel!,
            failure,input,slot;ndrange=length(input))
    end
    Backends.launch_kernel_async!(style,_compiled_paths_validate_amplitude_kernel!,
        failure,owner.inputs.pupil_amplitude,owner.reference;
        ndrange=length(owner.inputs.pupil_amplitude))
    return nothing
end

_compiled_paths_ratio!(::Backends.ScalarCPUStyle,attenuation,amplitude,reference,validation)=
    _compiled_paths_ratio!(attenuation,amplitude,reference)
KernelAbstractions.@kernel function _compiled_paths_ratio_kernel!(attenuation,amplitude,reference,failure)
    index=@index(Global,Linear)
    if index <= length(attenuation)
        a,r=amplitude[index],reference[index]
        value=zero(a)
        slot=0
        if !isfinite(a) || a < 0 || !isfinite(r) || r < 0
            slot=6
        elseif iszero(r)
            !iszero(a) && (slot=4)
        elseif iszero(a)
            slot=5
        else
            logratio=log(a)-log(r)
            if !isfinite(exp(logratio))
                slot=6
            else
                value=-logratio
            end
        end
        if !iszero(slot)
            KernelAbstractions.@atomic failure[slot]=UInt32(1)
        end
        attenuation[index]=value
    end
end
function _compiled_paths_ratio!(style::Backends.AcceleratorStyle,attenuation,amplitude,reference,validation)
    Backends.launch_kernel_async!(style,_compiled_paths_ratio_kernel!,
        attenuation,amplitude,reference,validation.device;ndrange=length(attenuation))
    return nothing
end

function _check_compiled_paths_inputs!(validation::CompiledPathsInputValidation)
    copyto!(parent(validation.host),validation.device)
    failure=validation.host
    iszero(failure[1]) || throw(ArgumentError("DM OPD must be finite"))
    iszero(failure[2]) || throw(ArgumentError("entrance OPD must be finite"))
    iszero(failure[3]) || throw(ArgumentError("entrance amplitude must be finite and nonnegative in reference-ratio mode"))
    iszero(failure[4]) || throw(ArgumentError("nonzero entrance amplitude where its reference is opaque"))
    iszero(failure[5]) || throw(ArgumentError("new opaque entrance samples inside reference support are unsupported; regenerate the optical model for changed pupil support"))
    iszero(failure[6]) || throw(ArgumentError("entrance amplitude ratio is not representable"))
    return nothing
end
_check_compiled_paths_inputs!(::Nothing)=nothing
_reset_compiled_paths_validation!(::Nothing)=nothing
function _reset_compiled_paths_validation!(validation::CompiledPathsInputValidation)
    fill!(validation.device,UInt32(0))
    fill!(validation.host,UInt32(0))
    return nothing
end
