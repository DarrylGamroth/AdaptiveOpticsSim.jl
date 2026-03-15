include(joinpath(@__DIR__, "common.jl"))

function transfer_functions(freq::AbstractVector{T}, loop_gain::T, Ti::T, Tau::T, Tdm::T) where {T<:AbstractFloat}
    S = complex.(zero(T), T(2π)) .* freq
    H_wfs = exp.(-Ti * S / 2)
    H_rtc = exp.(-Tau * S)
    H_dm = exp.(-Tdm * S)
    H_dac = (one(T) .- exp.(-S * Ti)) ./ (S * Ti)
    CC = loop_gain ./ (one(T) .- exp.(-S * Ti))
    H_ol = H_wfs .* H_rtc .* H_dac .* H_dm .* CC
    H_cl = H_ol ./ (one(T) .+ H_ol)
    H_er = one.(H_ol) ./ (one.(H_ol) .+ H_ol)
    H_n = H_cl ./ H_wfs
    return H_er, H_cl, H_ol, H_n
end

function main(; fs::Real=300.0, n_freq::Int=512, loop_gains=(0.2, 0.6))
    T = Float64
    freq = collect(range(T(fs / n_freq), T(fs / 2); length=n_freq - 1))
    Ti = T(1 / fs)
    Tau = Ti / 2
    Tdm = Ti / 2

    rejection_db = Matrix{T}(undef, length(freq), length(loop_gains))
    closed_loop_db = similar(rejection_db)
    for (idx, gain) in pairs(loop_gains)
        H_er, H_cl, _, _ = transfer_functions(freq, T(gain), Ti, Tau, Tdm)
        rejection_db[:, idx] .= 20 .* log10.(abs.(H_er))
        closed_loop_db[:, idx] .= 20 .* log10.(abs.(H_cl))
    end

    @info "Transfer-function tutorial complete" n_freq=length(freq) n_gains=length(loop_gains)
    return (
        frequency_hz=freq,
        rejection_db=rejection_db,
        closed_loop_db=closed_loop_db,
        loop_gains=collect(loop_gains),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
