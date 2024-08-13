using Gradus, Plots

# metric and metric parameters
m = KerrMetric(M=1.0, a=1.09)
# observer position
x = SVector(0.0, 1000.0, deg2rad(40), 0.0)
# accretion disc
d = ThinDisc(9.0, 60.0)

# define point function which filters geodesics that intersected the accretion disc
# and use those to calculate redshift
pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    2000.0,
    αlims = (-60, 60), 
    βlims = (-50, 50),
    image_width = 800,
    image_height = 400,
    callback = Gradus.domain_upper_hemisphere(),
    verbose = true,
    pf = pf,
)

heatmap(α, β, img, aspect_ratio = 1, colour = :magma, bgcolor =:black)