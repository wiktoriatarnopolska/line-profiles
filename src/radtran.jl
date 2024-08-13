using Gradus
using Plots
# metric and metric parameters

m = KerrMetric(M = 1.0, a = 0.5)
# observer position
x = SVector(0.0, 1000.0, deg2rad(80), 0.0)
# accretion disc
d = PolishDoughnut(m)

# set the emissivity
Gradus.emissivity_coefficient(::AbstractMetric, ::PolishDoughnut, x, ν) = 0.1

# define point function which reads the auxiliary variable
# which is contextually the intensity
pf = PointFunction((m, gp, t) -> gp.aux[1])

a, b, img = @time rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    2000.0,
    verbose = true,
    pf = pf,
    # instruct the integrator to solve the covariant radiative transfer equation
    αlims = (-25, 25), 
    βlims = (-15, 18),
    trace = Gradus.TraceRadiativeTransfer(I₀ = 0),
)

heatmap(a, b, img, aspect_ratio = 1, xlabel = "α", ylabel = "β", color =:gnuplot)
