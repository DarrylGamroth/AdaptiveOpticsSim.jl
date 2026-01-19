using Random

mutable struct Workspace{R<:AbstractRNG,T,A<:AbstractMatrix{Complex{T}},B<:AbstractMatrix{T}}
    rng::R
    pupil_field::A
    fft_buffer::A
    psf_buffer::B
end

function Workspace(n::Int; T=Float64, backend=Array, rng=MersenneTwister(0))
    pupil_field = backend{Complex{T}}(undef, n, n)
    fft_buffer = similar(pupil_field)
    psf_buffer = backend{T}(undef, n, n)
    return Workspace(rng, pupil_field, fft_buffer, psf_buffer)
end

function Workspace(ref::AbstractArray{S}, n::Int; T::Type{<:Real}=real(S), rng=MersenneTwister(0)) where {S}
    pupil_field = similar(ref, Complex{T}, n, n)
    fft_buffer = similar(pupil_field)
    psf_buffer = similar(ref, T, n, n)
    return Workspace(rng, pupil_field, fft_buffer, psf_buffer)
end

function ensure_psf_buffers!(ws::Workspace, n::Int)
    if size(ws.pupil_field, 1) != n || size(ws.pupil_field, 2) != n
        ws.pupil_field = similar(ws.pupil_field, n, n)
        ws.fft_buffer = similar(ws.pupil_field)
        ws.psf_buffer = similar(ws.psf_buffer, n, n)
    end
    return ws
end
