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

function _read_plane(file, index, ::Type{T}, shape) where {T}
    hdu = file[index]
    hdu isa ImageHDU{T,2} && size(hdu) == shape ||
        throw(ArgumentError("compiled optical array type or shape differs"))
    return read(hdu)
end

@inline _pixel(re, ::Nothing, i, j) = re[i,j]
@inline _pixel(re, im, i, j) = complex(re[i,j], im[i,j])

function _assemble_pixels(re::Matrix, im, ::Type{T}, swap) where {T}
    if !swap
        im === nothing && return eltype(re) === T ? re : T.(re)
        result = Matrix{T}(undef, size(re))
        result .= complex.(re, im)
        return result
    end
    result = Matrix{T}(undef, size(re, 2), size(re, 1))
    # A strip of input columns stays in cache while each output column is filled.
    for start in 1:32:size(re, 2)
        columns = start:min(start + 31, size(re, 2))
        for i in axes(re, 1), j in columns
            result[j,i] = _pixel(re, im, i, j)
        end
    end
    return result
end

function _read_pixels(file, real_hdu, imag_hdu, stored_type, output_type, shape, swap)
    real_hdu <= length(file) && imag_hdu <= length(file) ||
        throw(ArgumentError("compiled optical array HDU is missing"))
    disk_shape = swap ? reverse(shape) : shape
    re = _read_plane(file, real_hdu, real(stored_type), disk_shape)
    im = imag_hdu == 0 ? nothing : _read_plane(file, imag_hdu, real(stored_type), disk_shape)
    values = _assemble_pixels(re, im, output_type, swap)
    all(isfinite, values) || throw(ArgumentError(
        "compiled optical array contains nonfinite values or overflows conversion to $output_type"))
    return values
end

struct FITSArrayReader{F}
    decode::F
end

(reader::FITSArrayReader)(reference) = reader.decode(reference, nothing)
AdaptiveOpticsSim.Optics._read_compiled_array_as(reader::FITSArrayReader, reference, T) =
    reader.decode(reference, T)

function AOG._read_compiled_optical_model(path::AbstractString)
    path=realpath(path)
    yaml=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimYAMLExt)
    metadata=yaml._read_strict_yaml(path,ArgumentError)
    metadata isa AbstractDict && get(metadata,"schema_version",nothing)==1 || throw(ArgumentError("unsupported compiled optical model schema"))
    descriptors=get(metadata,"arrays",nothing)
    descriptors isa AbstractDict && !isempty(descriptors) || throw(ArgumentError("compiled optical model requires typed array descriptors"))
    root=dirname(path)
    function read_array(reference, output_type)
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
        return FITS(array_path, "r") do file
            _read_pixels(file, real_hdu, imag_hdu, T, something(output_type, T),
                (Int(shape[1]), Int(shape[2])), d["axes"][1] == 2)
        end
    end
    return metadata, FITSArrayReader(read_array)
end
end
