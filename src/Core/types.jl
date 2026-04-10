"""Base type for optical elements."""
abstract type AbstractOpticalElement end

"""Sources must provide wavelength(src)."""
abstract type AbstractSource <: AbstractOpticalElement end

"""Atmospheres implement advance!(atm, tel; rng) and propagate!(atm, tel)."""
abstract type AbstractAtmosphere <: AbstractOpticalElement end

"""
Wavefront sensors implement `measure!(wfs, tel[, src])` and
`update_valid_mask!(wfs, tel)`.

Optional behavior such as detector-coupled measurement, runtime preparation,
stacked-source support, and grouped execution is expressed through the
capability queries documented in the runtime/control interface layer.
"""
abstract type AbstractWFS <: AbstractOpticalElement end

function apply_shift_wfs!(::AbstractWFS; sx, sy)
    throw(InvalidConfiguration("apply_shift_wfs! is not supported for this WFS"))
end

"""Detectors implement capture!(det, psf; rng)."""
abstract type AbstractDetector <: AbstractOpticalElement end

"""Controllable optics implement command staging and `apply!(optic, tel, mode)`."""
abstract type AbstractControllableOptic <: AbstractOpticalElement end

abstract type DMApplyMode end
struct DMAdditive <: DMApplyMode end
struct DMReplace <: DMApplyMode end

struct RuntimeCommandSegment{L}
    label::L
    offset::Int
    length::Int
    function RuntimeCommandSegment(label, offset::Integer, length::Integer)
        offset >= 1 || throw(InvalidConfiguration("command segment offset must be >= 1"))
        length >= 0 || throw(InvalidConfiguration("command segment length must be >= 0"))
        return new{typeof(label)}(label, Int(offset), Int(length))
    end
end

struct RuntimeCommandLayout{S}
    segments::S
    total_length::Int
    function RuntimeCommandLayout(segments::Tuple)
        total = 0
        for seg in segments
            seg isa RuntimeCommandSegment ||
                throw(InvalidConfiguration("RuntimeCommandLayout segments must be RuntimeCommandSegment values"))
            seg.offset == total + 1 ||
                throw(InvalidConfiguration("RuntimeCommandLayout segments must be contiguous and 1-based"))
            total += seg.length
        end
        return new{typeof(segments)}(segments, total)
    end
end

RuntimeCommandLayout(segments::RuntimeCommandSegment...) = RuntimeCommandLayout(tuple(segments...))
RuntimeCommandLayout(specs::Pair...) = begin
    offset = 1
    segments = map(specs) do spec
        length = Int(spec.second)
        seg = RuntimeCommandSegment(spec.first, offset, length)
        offset += length
        seg
    end
    RuntimeCommandLayout(Tuple(segments))
end
command_segments(layout::RuntimeCommandLayout) = layout.segments
command_segment_labels(layout::RuntimeCommandLayout) = map(seg -> seg.label, layout.segments)
command_segment_range(seg::RuntimeCommandSegment) = seg.offset:(seg.offset + seg.length - 1)

"""Deformable mirrors implement build_influence_functions!(dm, tel) and apply!(dm, tel, mode)."""
abstract type AbstractDeformableMirror <: AbstractControllableOptic end
