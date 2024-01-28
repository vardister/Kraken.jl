"""
    wilson_charge_jb(z, fs, T)

Generate a SUS charge signature according to Wilson et al. 2021:
https://ieeexplore.ieee.org/document/8835907/

Required arguments:
- `z`: source depth (m)
- `fs`: sampling frequency (Hz)
- `T`: total duration (s)

Return:
- `Pfinal`: pressure signature (MPa)

"""
function wilson_charge_jb(z, fs, T)
    R = 1.0
    z_0 = z + 10.1
    w = 0.03118 * 1.25
    w_third = w^(1/3)
    wR_ratio = w_third/R

    P_s = 5.04e13 * (wR_ratio)^(1.13)
    tau_s = 8.12e-5 * w_third * (wR_ratio)^(- .14)
    T1 = 2 * w_third * z_0^(-5/6)
    T2 = .8 * T1
    theta_s = 0.194 * w_third * z_0^(-5/6)
    P_1 = 1.49e12 * wR_ratio * z_0^(.33)
    P_2 = 3.93e11 * wR_ratio * z_0^(.28)
    P_min1 = 5e10 * wR_ratio * z_0^(.6)
    P_min2 = .58 * P_min1
    tau_r = 1.36e-2 * w_third * z_0^(-0.6)
    tau_d = 0.87e-2 * w_third * z_0^(-0.6)
    tau_r2 = 2.5 * tau_r
    tau_d2 = 2.5 * tau_d
    theta_1 = 0.45 * w_third * z_0^(-5/6)
    theta_min1 = 1.64 * w_third * z_0^(-5/6)
    theta_min2 = .69 * theta_min1

    dt = 1 / fs
    t = 0:dt:0.25
    t_max = maximum(t)

    Pa = P_s .* exp.(-t ./ tau_s)
    Pa = [0; Pa]
    Pa = Pa[1:end-1]

    Pa_short = copy(Pa)
    max_Pa = maximum(Pa)
    Pa_short[Pa_short .< 0.1 * max_Pa] .= 0
    xfind = findlast(x -> x > 0.1 * max_Pa, Pa)
    t_short = t[xfind]
    Pa_long = 0.1 * max_Pa .* exp.(-(t .- t_short) ./ (15 * tau_s))
    Pa_long[t .<= t_short] .= 0
    Pa_2scale = Pa_short + Pa_long

    xx = findfirst(t -> t > theta_s + dt, t)
    fnegphase_1 = 1 / (2 * theta_min1)
    Pb1 = -P_min1 .* sin.(2 * π * fnegphase_1 .* (t .- t[xx]))
    Pb1[Pb1 .> 0] .= 0
    Pb1[t .+ dt .> theta_s + theta_min1] .= 0

    Pc1 = P_1 .* exp.((t .- T1) ./ tau_r)
    Pc1[t .> T1 + dt] .= 0
    Pd1 = P_1 .* exp.(-(t .- T1) ./ tau_d)
    Pd1[t .< T1 + dt] .= 0

    P_1_bubblepulse = Pa_2scale + Pb1 + Pc1 + Pd1

    Pc2 = P_2 .* exp.((t .- (T1 + T2)) ./ tau_r2)
    Pc2[Pc2 .> P_2] .= 0

    Pd2 = P_2 .* exp.(-(t .- (T1 + T2)) ./ tau_d2)
    Pd2[t .< (T1 + T2)] .= 0

    xfind = findfirst(t -> t >= dt + theta_s + theta_min1 + theta_1, t)
    fnegphase_2 = 1 / (2 * theta_min2)
    Pb2 = -P_min2 .* sin.(2 * π * fnegphase_2 .* (t .- t[xfind]))
    Pb2[t .< dt + theta_s + theta_min1 + theta_1] .= 0
    Pb2[t .> dt + theta_s + theta_min1 + theta_1 + theta_min2] .= 0

    Pfinal = P_1_bubblepulse + Pb2 + Pc2 + Pd2

    tminus = -500:-1
    tminus *= dt
    Pminus = zeros(length(tminus))

    t = [tminus; t]
    P_1_bubblepulse = [Pminus; P_1_bubblepulse]
    Pfinal = [Pminus; Pfinal]

    append!(Pfinal, zeros(Float64, T*fs - length(Pfinal)))
    Pfinal /= 1e6

    return Pfinal
end
