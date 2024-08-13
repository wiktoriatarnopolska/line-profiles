using Gradus, Plots

m1 = KerrMetric(M = 1.0, a = 0.998)
x = SVector(0.0, 10000.0, deg2rad(40), 0.0)
d = ThinDisc(9.0, 60.0)


α, β, img1 = rendergeodesics(
    m1,
    x,
    d,
    # max integration time
    20_000.0,
    image_width = 800,
    image_height = 800,
    αlims = (-30, 30),
    βlims = (-20, 20),
    verbose = true,
    ensemble = Gradus.EnsembleEndpointThreads(),
    callback = Gradus.domain_upper_hemisphere(),
    pf = ConstPointFunctions.redshift(m1, x) ∘ ConstPointFunctions.filter_intersected(),
)

m2 = KerrMetric(M = 1.0, a = -0.998)

α, β, img2 = rendergeodesics(
    m2,
    x,
    d,
    # max integration time
    20_000.0,
    image_width = 800,
    image_height = 800,
    αlims = (-30, 30),
    βlims = (-20, 20),
    verbose = true,
    ensemble = Gradus.EnsembleEndpointThreads(),
    callback = Gradus.domain_upper_hemisphere(),
    pf = ConstPointFunctions.redshift(m2, x) ∘ ConstPointFunctions.filter_intersected(),
)

p1 = heatmap(
    α, 
    β,
    img1,
    color = :magma,
    xlabel = "α",
    ylabel = "β",
    aspect_ratio = 1,
    minorgrid = true,
    title = "a = 0.998",
    bgcolour = :black,

)
# contour!(p1, α, β, img1, color = :red)

p2 = heatmap(
    α, 
    β,
    img2,
    color = :magma,
    xlabel = "α",
    ylabel = "β",
    aspect_ratio = 1,
    minorgrid = true,
    title = "a = -0.998",
    bgcolour = :black,

)
# contour!(p2, α, β, img2, color = :red)

p = heatmap(
    α,
    β,
    img1 .-img2,
    color = :magma,
    xlabel = "α",
    ylabel = "β",
    aspect_ratio = 1,
    minorgrid = true,
    title = "energy shift",
    bgcolour = :black,
)
# contour!(α, β, img2 .- img1, color = :red)
