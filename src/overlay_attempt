using Gradus
using Plots
using Printf
using Images
gr()

_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.3f" m.a

_format_label(edd) = Printf.@sprintf "Ṁ / Ṁedd = %.1f" (edd / 100)

function calculate_line_profile(m, x, d, prof, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d,
        prof;
        method = TransferFunctionMethod(),
        # try binning???
        minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        #maxrₑ = 500.0,
        maxrₑ = 30.0,
        # resolution
        numrₑ = 200,
        #numrₑ = 50,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end

function run_all_parameter_combinations(model, θ, bins; kwargs...)

    x = SVector(0.0, 1000.0, deg2rad(θ), 0.0)

    # You can change the spin here
    m = KerrMetric(M = 1.0, a = 0.99)

    spectrum = PowerLawSpectrum(2.0)

    @info "m = $m θ = $(θ)"

    # discs
    dthin = ThinDisc(Gradus.isco(m), Inf)
    d1 = ShakuraSunyaev(m, eddington_ratio = 0.1)
    d2 = ShakuraSunyaev(m, eddington_ratio = 0.2)
    d3 = ShakuraSunyaev(m, eddington_ratio = 0.3)

    @info "0%"
    em_prof0 = em_prof(m, dthin, model, spectrum)
    prof0 = Gradus.RadialDiscProfile(em_prof0)
    edrat0 = @time calculate_line_profile(m, x, dthin, prof0, bins; kwargs...)
    @info "10%"
    em_prof10 = em_prof(m, d1, model, spectrum)
    prof10 = Gradus.RadialDiscProfile(em_prof10)
    edrat10 = @time calculate_line_profile(m, x, d1, prof10, bins; kwargs...)
    @info "20%"
    em_prof20 = em_prof(m, d2, model, spectrum)
    prof20 = Gradus.RadialDiscProfile(em_prof20)
    edrat20 = @time calculate_line_profile(m, x, d2, prof20, bins; kwargs...)
    @info "30%"
    em_prof30 = em_prof(m, d3, model, spectrum)
    prof30 = Gradus.RadialDiscProfile(em_prof30)
    edrat30 = @time calculate_line_profile(m, x, d3, prof30, bins; kwargs...)

    return (; metric = m, f = [edrat0, edrat10, edrat20, edrat30], θ = θ, bins = bins, h = model.h)

end


bins = collect(range(0.1, 1.4, 200))

function em_prof(m, d, model, spectrum; kwargs...)
    em_prof = Gradus.emissivity_profile(
        m,
        d,
        model,
        spectrum,
        n_samples = 10_000,
    )
    return em_prof
end

# models
model3 = LampPostModel(h = 3.0)
model6 = LampPostModel(h = 6.0)
model12 = LampPostModel(h = 12.0)

# inclinations
inclination15 = 15.0
inclination30 = 30.0
inclination60 = 60.0

#moving offset 0 -> 1 to avoid transfer functions' confusion
KWARGS = (; β₀ = 1)

# inclination = 15
data_n1 = run_all_parameter_combinations(model3, inclination15, bins; KWARGS...)

function plot_all(data, p)
    #  incl_text = Printf.@sprintf " θ=%0.f h=%0.f" data.θ data.h
    # p = plot(title = _format_metric(data.metric) * incl_text, legend = :topleft)
    m = Float64[1.0, 1.0, 1.0, 1.0]

    m[1] = maximum(data.f[1])
    m[2] = 1.06*m[1]
    m[3] = 1.16*m[1]
    m[4] = 1.35*m[1]
    for (edd, f, i) in zip((0, 10, 20, 30), data.f, (1, 2, 3, 4))
        plot!(p, data.bins, f/m[i], label = _format_label(edd), linewidth = 3)
    end
    p
end

img = load("test_plot.png")
(x_scale, y_scale, ar) = (-0.34432783608195905:0.003711550817997595:1.6042363433667783, -0.279629126185809:0.002474967061668832:1.0939775930403928, 1.0570342205323193)
p = plot(x_scale, y_scale, reverse(img, dims = 1), yflip = false, aspect_ratio = ar, size=(800,800))

plot_all(data_n1, p)


# inclination = 30

data_n1 = run_all_parameter_combinations(model3, inclination30, bins; KWARGS...)

function plot_all(data, p)
    #  incl_text = Printf.@sprintf " θ=%0.f h=%0.f" data.θ data.h
    # p = plot(title = _format_metric(data.metric) * incl_text, legend = :topleft)
    m = Float64[1.0, 1.0, 1.0, 1.0]

    m[1] = maximum(data.f[1])
    m[2] = 1.06*m[1]
    m[3] = 1.16*m[1]
    m[4] = 1.35*m[1]
    for (edd, f, i) in zip((0, 10, 20, 30), data.f, (1, 2, 3, 4))
        plot!(p, data.bins, f/m[i], label = _format_label(edd), linewidth = 3)
    end
    p
end

img = load("test_plot2.png")
(x_scale, y_scale, ar) = (-0.34432783608195905:0.003711550817997595:1.6042363433667783, -0.279629126185809:0.002474967061668832:1.0939775930403928, 1.0570342205323193)
p = plot(x_scale, y_scale, reverse(img, dims = 1), yflip = false, aspect_ratio = ar, size=(800,800))

plot_all(data_n1, p)


# inclination = 60

data_n1 = run_all_parameter_combinations(model3, inclination60, bins; KWARGS...)

function plot_all(data, p)
    #  incl_text = Printf.@sprintf " θ=%0.f h=%0.f" data.θ data.h
    # p = plot(title = _format_metric(data.metric) * incl_text, legend = :topleft)
    m = Float64[1.0, 1.0, 1.0, 1.0]

    m[1] = maximum(data.f[1])
    m[2] = 1.06*m[1]
    m[3] = 1.16*m[1]
    m[4] = 1.35*m[1]
    for (edd, f, i) in zip((0, 10, 20, 30), data.f, (1, 2, 3, 4))
        plot!(p, data.bins, f/m[i], label = _format_label(edd), linewidth = 3)
    end
    p
end

img = load("test_plot3.png")
#(x_scale, y_scale, ar) = (-0.34432783608195905:0.003711550817997595:1.6042363433667783, -0.279629126185809:0.002474967061668832:1.0939775930403928, 1.0570342205323193)
p = plot(x_scale, y_scale, reverse(img, dims = 1), yflip = false, aspect_ratio = ar, size=(800,800))

plot_all(data_n1, p)