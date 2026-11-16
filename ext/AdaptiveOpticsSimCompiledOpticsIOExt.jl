module AdaptiveOpticsSimCompiledOpticsIOExt
using AdaptiveOpticsSim
using FITSIO
using YAML
import AdaptiveOpticsSim.AlgorithmGraphs as AOG

function _asset(root,relative)
    relative isa String && !isempty(relative) && !isabspath(relative) && normpath(relative)==relative ||
        throw(ArgumentError("compiled optical asset requires a normalized relative path"))
    all(x->!(x in (".","..")),splitpath(relative)) || throw(ArgumentError("compiled optical asset escapes model directory"))
    path=realpath(joinpath(root,relative))
    startswith(path,root*"/") || throw(ArgumentError("compiled optical asset escapes model directory"))
    return path
end
function AOG._read_compiled_optical_model(path::AbstractString)
    path=realpath(path)
    yaml=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimYAMLExt)
    metadata=yaml._read_strict_yaml(path,ArgumentError)
    metadata isa AbstractDict && get(metadata,"schema_version",nothing)==1 || throw(ArgumentError("unsupported compiled optical model schema"))
    descriptors=get(metadata,"arrays",nothing)
    descriptors isa AbstractDict && !isempty(descriptors) || throw(ArgumentError("compiled optical model requires typed array descriptors"))
    root=dirname(path)
    function read_array(reference)
        haskey(descriptors,reference) || throw(ArgumentError("compiled optical array reference is undeclared"))
        d=descriptors[reference]
        required = ("path","shape","type","axes","real_hdu","imag_hdu")
        all(key -> haskey(d, key), required) &&
            all(key -> key in required || key == "sha256", keys(d)) ||
            throw(ArgumentError("unknown or missing compiled optical array descriptor fields"))
        array_path=_asset(root,d["path"])
        types=Dict("Float32"=>Float32,"Float64"=>Float64,"ComplexF32"=>ComplexF32,"ComplexF64"=>ComplexF64)
        haskey(types,d["type"]) || throw(ArgumentError("unsupported compiled optical array type"))
        d["axes"] isa AbstractVector && d["axes"] in ([1,2],[2,1]) &&
            all(x->x isa Integer && !(x isa Bool),d["axes"]) || throw(ArgumentError("compiled optical array axes must permute[y,x]"))
        shape=d["shape"]
        shape isa AbstractVector && length(shape)==2 && all(x->x isa Integer && !(x isa Bool) && x>0,shape) ||
            throw(ArgumentError("compiled optical arrays require a positive rank-two shape"))
        real_hdu,imag_hdu=d["real_hdu"],d["imag_hdu"]
        real_hdu isa Integer && !(real_hdu isa Bool) && real_hdu>=1 || throw(ArgumentError("invalid real_hdu"))
        imag_hdu isa Integer && !(imag_hdu isa Bool) && imag_hdu>=0 || throw(ArgumentError("invalid imag_hdu"))
        T=types[d["type"]]
        (T<:Complex)==(imag_hdu>0) || throw(ArgumentError("compiled optical array type disagrees with imaginary HDU"))
        values=FITS(array_path,"r") do file
            real_hdu<=length(file) && imag_hdu<=length(file) || throw(ArgumentError("compiled optical array HDU is missing"))
            re=permutedims(read(file[real_hdu]),Tuple(d["axes"]))
            imag_hdu==0 ? re : complex.(re,permutedims(read(file[imag_hdu]),Tuple(d["axes"])))
        end
        collect(size(values))==shape && eltype(values)===T && all(isfinite,values) ||
            throw(ArgumentError("compiled optical array type, shape or finite-value contract differs"))
        return values
    end
    return metadata,read_array
end
end
