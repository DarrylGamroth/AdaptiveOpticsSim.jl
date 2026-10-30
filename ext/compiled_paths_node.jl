# Generic artifact adapter for the lowered multi-path Fresnel executor. Its live
# input grids are independent of propagation and fitted-screen sample counts.
struct CompiledPathsNode{T} end
struct CompiledPathsConfiguration{D,P,R,T}
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
    c=CompiledPathsConfiguration(path,hash,dm_grid,pupil_grid,reference,schemas,mode,amplitude_mode,
        summation,Int(metadata["image_size"]),T(metadata["dx_output"]))
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
struct CompiledPathsInputWorkspace{A}
    amplitude::A
    reference::A
end
struct CompiledPathsOwner{P,I,O,M,W,R,T}
    propagation::P
    inputs::I
    outputs::O
    mapping::M
    workspace::W
    reference::R
    sampling_m::T
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
    metadata,read_array=AOG._read_compiled_optical_model(c.path,c.sha256)
    live=metadata["live"]
    ndm=Int(live["dm_size"])
    nentrance=Int(live["entrance_size"])
    dm=zeros(T,ndm,ndm)
    entrance=zeros(T,nentrance,nentrance)
    attenuation=similar(entrance);fill!(attenuation,0)
    propagation=_prepare_compiled_paths(metadata,read_array,nothing;T,target,dm_opd=dm,
        entrance_opd=entrance,entrance_attenuation=attenuation,summation=c.summation)
    dm_mapping=_bilinear_interpolation_plan(c.dm_grid,ndm,T(live["dm_pitch"]),zero(T))
    pupil_mapping=c.mode===:sky ? _bilinear_interpolation_plan(c.pupil_grid,nentrance,T(live["entrance_pitch"]),zero(T)) : nothing
    native_reference=similar(entrance)
    native_amplitude=similar(entrance)
    if c.mode===:sky && c.reference!==nothing
        _interpolate_grid!(native_reference,c.reference,pupil_mapping)
    else
        fill!(native_reference,1)
    end
    owner=CompiledPathsOwner(propagation,inputs,outputs,(dm=dm_mapping,pupil=pupil_mapping),
        CompiledPathsInputWorkspace(native_amplitude,native_reference),c.reference,c.sampling_m)
    return owner
end
function AOG.step_graph_node!(owner::CompiledPathsOwner)
    p=owner.propagation
    _compiled_paths_validate_inputs!(owner,owner.reference)
    _interpolate_grid!(p.inputs.dm_opd,owner.inputs.dm_opd,owner.mapping.dm)
    _compiled_paths_entrance!(owner,owner.mapping.pupil)
    products=_enqueue_compiled_paths!(p)
    for name in (:fringed,:unfringed,:reference)
        copyto!(getproperty(owner.outputs,name),getproperty(products,name))
    end
    return nothing
end
_compiled_paths_entrance!(owner,::Nothing)=nothing
function _compiled_paths_entrance!(owner,mapping::BilinearInterpolationPlan)
    inputs,p=owner.inputs,owner.propagation
    all(isfinite,inputs.pupil_opd) || throw(ArgumentError("entrance OPD must be finite"))
    _interpolate_grid!(p.inputs.entrance_opd,inputs.pupil_opd,mapping)
    _compiled_paths_amplitude!(owner,owner.reference,mapping)
    return nothing
end
function _compiled_paths_amplitude!(owner,::Nothing,mapping)
    a=owner.inputs.pupil_amplitude
    all(isfinite,a) || throw(ArgumentError("entrance log attenuation must be finite"))
    _interpolate_grid!(owner.propagation.inputs.entrance_attenuation,a,mapping)
    return nothing
end
function _compiled_paths_amplitude!(owner,reference::AbstractMatrix,mapping)
    a=owner.inputs.pupil_amplitude
    _interpolate_grid!(owner.workspace.amplitude,a,mapping)
    _compiled_paths_ratio!(owner.propagation.inputs.entrance_attenuation,owner.workspace.amplitude,owner.workspace.reference)
    return nothing
end
function AOG.reset_graph_node!(owner::CompiledPathsOwner)
    for name in (:fringed,:unfringed,:reference)
        fill!(getproperty(owner.outputs,name),0)
    end
    return nothing
end
AOG.optical_sampling_m(owner::CompiledPathsOwner)=owner.sampling_m
AOG.graph_node_capture_capability(owner::CompiledPathsOwner)=AOG.GraphNodeCaptureUnsupported()

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
_require_compiled_input_target(target)=throw(AOG.AlgorithmGraphError("compiled optical input-grid adapter currently supports CPU execution"))
