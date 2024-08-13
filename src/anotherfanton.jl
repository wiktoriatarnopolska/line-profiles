using Gradus
using Plots

# define custom bins for g
bins = collect(range(0.1, 2.0, 300))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)

# inclination
x = SVector(0.0, 1000.0, deg2rad(30), 0.0)


#geometric thin disc
d = ThinDisc(0.0, Inf)

m = KerrMetric(M = 1.0, a = 0.998)

function calculate_line_profile1(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        minrₑ = 1.23,
        maxrₑ = 10.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end

function calculate_line_profile2(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        minrₑ = 1.23,
        maxrₑ = 20.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end

function calculate_line_profile3(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        minrₑ = 1.23,
        maxrₑ = 50.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end

function calculate_line_profile4(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        minrₑ = 1.23,
        maxrₑ = 100.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end

edrat1 = calculate_line_profile1(m, x, d, bins)
edrat2 = calculate_line_profile2(m, x, d, bins)
edrat3 = calculate_line_profile3(m, x, d, bins)
edrat4 = calculate_line_profile4(m, x, d, bins)

plot(bins, edrat1, label = "r_out = 10", title = "a = 0.998, r_in = 1.23")
plot!(bins, edrat2, label = "r_out = 20")
plot!(bins, edrat3, label = "r_out = 50")
plot!(bins, edrat4, label = "r_out = 100")
