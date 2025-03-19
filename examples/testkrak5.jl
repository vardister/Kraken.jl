using Kraken
using CairoMakie

frq=3:100
nm = 20
nl = 1
note1 = b"NVW"
note2 = b"A"

b1 = [0.0 0.0 100.0]
ssp1 = [0.000 1500.0 0.000 1.030 0.000 0.000
        100.0 1500.0 0.000 1.030 0.000 0.000]

bsig = 0.0
clh = [0.0 1800.0]
rng = 5000.0
nsr = 1
zsr = 1.0
nrc = 101
zrc = [0.0, 100.0]

b = b1
ssp = ssp1
nc, iq = size(ssp)

sspTHS = [0.0   343.0 0.0 0.00121 0.00 0.0]
sspBHS = [100.0 1800.0 0.0 2.0 0.00 0.0]
sspHS = [sspTHS; sspBHS]

nz = nsr + nrc

env = EnvKRAKEN(
    ssp,
    b,
    sspHS,
    zrc,
    zsr,
)

#%% Run KRAKEN
Nz = nrc + 1 # the number of depth points
CP = zeros(nm, length(frq))
CG = zeros(nm, length(frq))
kr_real = zeros(nm, length(frq))
kr_imag = zeros(nm, length(frq))
phi = zeros(Nz, nm, length(frq))

for (i, ff) in enumerate(frq)
    res = kraken(env, ff; n_modes = nm)
    CP[:, i] = res["cp"]
    CG[:, i] = res["cg"]
    kr_real[:, i] = res["kr_real"]
    kr_imag[:, i] = res["kr_imag"]
    phi[:, :, i] = res["modes"]
end

#%% Plot group and phase speeds
# before plotting, replace all zeros with NaN
CP[CP .== 0] .= NaN
CG[CG .== 0] .= NaN

f = Figure();
ax1 = Axis(f[1, 1], xlabel = "Frequency (m/s)", ylabel = "Phase Speed (m/s)")
series!(ax1, frq, CP; color=:tab20)
ax2 = Axis(f[1, 2], ylabel = "Group Speed (m/s)", xlabel="Frequency (Hz)")
series!(ax2, frq, CG; color=:tab20)
display(f)

# #%% Plot modes
# # before plotting, replace all zeros with NaN
# phi[phi .== 0] .= NaN
# f2 = Figure()
# ax1 = Axis(f2[1, 1], ylabel = "Depth (m)", xlabel = "Mode Amplitude")
# for i in 1:nm
#     lines!(ax1, phi[:, i, 39], -vec(res["zm"]))
# end
# display(f2)

#%% Plot horizontal wavenumbers
# # before plotting, replace all zeros with NaN
kr_real[kr_real .== 0] .= NaN
kr_imag[kr_imag .== 0] .= NaN

f3 = Figure()
ax1 = Axis(f3[1, 1], ylabel = "Frequency (Hz)", xlabel = "Wavenumber - real part",
           xticklabelrotation = pi/4)
for i in 1:nm
    lines!(ax1, kr_real[i, :], frq)
end

ax2 = Axis(f3[1, 2], xlabel = "Wavenumber - imaginary part",
           xticklabelrotation = pi/4)
for i in 1:nm
    lines!(ax2, kr_imag[i, :], frq)
end
hideydecorations!(ax2, grid = false)
display(f3)