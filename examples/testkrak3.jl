using KRAKEN
using CairoMakie

frq = 150.0
nm = 10
nl = 2
note1 = b"SVW"
b1 = [9000.0 0.0 5000.000]

ssp1 = [
    0.000 1509.520 0.000 1.03000 0.000 0.000
    10.000 1509.500 0.000 1.03000 0.000 0.000
    20.000 1509.490 0.000 1.03000 0.000 0.000
    30.000 1509.290 0.000 1.03000 0.000 0.000
    40.000 1508.520 0.000 1.03000 0.000 0.000
    50.000 1507.240 0.000 1.03000 0.000 0.000
    60.000 1505.640 0.000 1.03000 0.000 0.000
    70.000 1503.920 0.000 1.03000 0.000 0.000
    80.000 1502.290 0.000 1.03000 0.000 0.000
    90.000 1500.810 0.000 1.03000 0.000 0.000
    100.000 1499.390 0.000 1.03000 0.000 0.000
    110.000 1497.930 0.000 1.03000 0.000 0.000
    120.000 1496.460 0.000 1.03000 0.000 0.000
    130.000 1495.020 0.000 1.03000 0.000 0.000
    140.000 1493.670 0.000 1.03000 0.000 0.000
    150.000 1492.420 0.000 1.03000 0.000 0.000
    160.000 1491.310 0.000 1.03000 0.000 0.000
    170.000 1490.330 0.000 1.03000 0.000 0.000
    180.000 1489.480 0.000 1.03000 0.000 0.000
    190.000 1488.740 0.000 1.03000 0.000 0.000
    200.000 1488.110 0.000 1.03000 0.000 0.000
    210.000 1487.580 0.000 1.03000 0.000 0.000
    220.000 1487.120 0.000 1.03000 0.000 0.000
    230.000 1486.710 0.000 1.03000 0.000 0.000
    240.000 1486.340 0.000 1.03000 0.000 0.000
    250.000 1485.980 0.000 1.03000 0.000 0.000
    260.000 1485.620 0.000 1.03000 0.000 0.000
    270.000 1485.250 0.000 1.03000 0.000 0.000
    280.000 1484.890 0.000 1.03000 0.000 0.000
    290.000 1484.540 0.000 1.03000 0.000 0.000
    300.000 1484.210 0.000 1.03000 0.000 0.000
    310.000 1483.890 0.000 1.03000 0.000 0.000
    320.000 1483.610 0.000 1.03000 0.000 0.000
    330.000 1483.340 0.000 1.03000 0.000 0.000
    340.000 1483.090 0.000 1.03000 0.000 0.000
    350.000 1482.860 0.000 1.03000 0.000 0.000
    360.000 1482.640 0.000 1.03000 0.000 0.000
    375.000 1482.340 0.000 1.03000 0.000 0.000
    400.000 1481.890 0.000 1.03000 0.000 0.000
    425.000 1481.490 0.000 1.03000 0.000 0.000
    450.000 1481.140 0.000 1.03000 0.000 0.000
    475.000 1480.840 0.000 1.03000 0.000 0.000
    500.000 1480.600 0.000 1.03000 0.000 0.000
    525.000 1480.430 0.000 1.03000 0.000 0.000
    550.000 1480.310 0.000 1.03000 0.000 0.000
    575.000 1480.230 0.000 1.03000 0.000 0.000
    600.000 1480.200 0.000 1.03000 0.000 0.000
    650.000 1480.200 0.000 1.03000 0.000 0.000
    700.000 1480.270 0.000 1.03000 0.000 0.000
    750.000 1480.400 0.000 1.03000 0.000 0.000
    800.000 1480.570 0.000 1.03000 0.000 0.000
    850.000 1480.770 0.000 1.03000 0.000 0.000
    900.000 1481.000 0.000 1.03000 0.000 0.000
    950.000 1481.250 0.000 1.03000 0.000 0.000
    1000.000 1481.530 0.000 1.03000 0.000 0.000
    1050.000 1481.830 0.000 1.03000 0.000 0.000
    1100.000 1482.170 0.000 1.03000 0.000 0.000
    1150.000 1482.530 0.000 1.03000 0.000 0.000
    1200.000 1482.900 0.000 1.03000 0.000 0.000
    1250.000 1483.270 0.000 1.03000 0.000 0.000
    1300.000 1483.640 0.000 1.03000 0.000 0.000
    1350.000 1484.040 0.000 1.03000 0.000 0.000
    1400.000 1484.490 0.000 1.03000 0.000 0.000
    1450.000 1484.970 0.000 1.03000 0.000 0.000
    1500.000 1485.470 0.000 1.03000 0.000 0.000
    1550.000 1485.960 0.000 1.03000 0.000 0.000
    1600.000 1486.450 0.000 1.03000 0.000 0.000
    1650.000 1486.930 0.000 1.03000 0.000 0.000
    1700.000 1487.430 0.000 1.03000 0.000 0.000
    1750.000 1487.950 0.000 1.03000 0.000 0.000
    1800.000 1488.510 0.000 1.03000 0.000 0.000
    1850.000 1489.090 0.000 1.03000 0.000 0.000
    1900.000 1489.700 0.000 1.03000 0.000 0.000
    1950.000 1490.340 0.000 1.03000 0.000 0.000
    2000.000 1490.990 0.000 1.03000 0.000 0.000
    2100.000 1492.350 0.000 1.03000 0.000 0.000
    2200.000 1493.760 0.000 1.03000 0.000 0.000
    2300.000 1495.220 0.000 1.03000 0.000 0.000
    2400.000 1496.710 0.000 1.03000 0.000 0.000
    2500.000 1498.230 0.000 1.03000 0.000 0.000
    2750.000 1502.140 0.000 1.03000 0.000 0.000
    3000.000 1506.180 0.000 1.03000 0.000 0.000
    3250.000 1510.330 0.000 1.03000 0.000 0.000
    3500.000 1514.530 0.000 1.03000 0.000 0.000
    3750.000 1518.760 0.000 1.03000 0.000 0.000
    4000.000 1523.210 0.000 1.03000 0.000 0.000
    4250.000 1527.840 0.000 1.03000 0.000 0.000
    4500.000 1531.450 0.000 1.03000 0.000 0.000
    4750.000 1533.370 0.000 1.03000 0.000 0.000
    5000.000 1536.150 0.000 1.03000 0.000 0.000
]

b2 = [600.0 0.0 5300.0]
ssp2 = [
    5000.000 1513.108 0.0 1.75100 0.000 0.000
    5300.000 1613.108 0.0 1.75100 0.000 0.000
]
note2 = b"R"
bsig = 0.0
clh = [1300.0 1540.0]
rng = 5000000.0
nsr = 1
zsr = 500.0
nrc = 500
zrc = [0.0 5300.0]

b = [b1; b2]
ssp = [ssp1; ssp2]
nc, iq = size(ssp)

sspTHS = [0.0 1700.0 0.0 1.7 0.0 0.0]
sspBHS = [5300.0 1749.0 0.0 1.941 0.1 0.0]
sspHS = [sspTHS; sspBHS]

nz = nsr + nrc
frq=vcat(10:1:30, 32:2:50, 55:5:100, 110:10:250)

env = Env(
    ssp = ssp,
    b = b,
    sspHS = sspHS,
    n_krak = nrc,
    z_krak = zrc,
    z_sr = zsr,
    n_sr = nsr,
)

CP=zeros(nm, length(frq))
CG=zeros(nm, length(frq))

for (i, ff) in enumerate(frq)
    res = kraken(env, ff; n_modes = nm)
    CP[:, i] = res["cp"]
    CG[:, i] = res["cg"]
end

#%% Plot group and phase speeds
f = Figure()
ax1 = Axis(f[1, 1], xlabel = "Frequency (Hz)", ylabel = "Phase Speed (m/s)", title = "Phase Speed")
series!(ax1, frq, CP, color = :tab20)
ax2 = Axis(f[1, 2], xlabel = "Frequency (Hz)", ylabel = "Group Speed (m/s)", title = "Group Speed")
series!(ax2, frq, CG, color = :tab20)
display(f)
