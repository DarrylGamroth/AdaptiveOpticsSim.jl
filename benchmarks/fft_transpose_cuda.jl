# Benchmark-only coalesced transpose; no production FFT policy depends on it.
module FFTTransposeExperiment
import CUDA
using KernelAbstractions
const TRANSPOSE = get(ENV, "AOS_FFT_TRANSPOSE", "generic")
TRANSPOSE in ("generic", "tiled") || error("AOS_FFT_TRANSPOSE must be generic or tiled")

@kernel function tiled_transpose_kernel!(destination, source)
    x, y, _ = @index(Local, NTuple)
    bx, by, plane = @index(Group, NTuple)
    tile = @localmem eltype(source) (33, 32)
    source_x = (bx-1)*32+x
    source_y = (by-1)*32+y
    for offset in 0:8:24
        if source_x <= size(source,1) && source_y+offset <= size(source,2)
            @inbounds tile[x,y+offset] = source[source_x,source_y+offset,plane]
        end
    end
    @synchronize
    destination_x = (by-1)*32+x
    destination_y = (bx-1)*32+y
    for offset in 0:8:24
        if destination_x <= size(destination,1) && destination_y+offset <= size(destination,2)
            @inbounds destination[destination_x,destination_y+offset,plane] = tile[y+offset,x]
        end
    end
end

function transpose_planes!(destination, source)
    if TRANSPOSE == "generic"
        return permutedims!(destination, source, (2,1,3))
    end
    kernel = tiled_transpose_kernel!(get_backend(source), (32,8,1))
    kernel(destination, source; ndrange=(cld(size(source,1),32)*32,
        cld(size(source,2),32)*8, size(source,3)))
    return destination
end
end
