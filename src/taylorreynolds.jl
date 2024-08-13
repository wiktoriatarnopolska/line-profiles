using Gradus
using Plots
using Printf
using LaTeXStrings
using Plots.PlotMeasures

_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.3f" m.a

_format_label(edd) = Printf.@sprintf "Ṁ / Ṁedd = %.1f" (edd / 100)

function calculate_line_profile(m, x, d, prof, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d,
        prof;
        method = BinningMethod(),
        # try binning???
        minrₑ = Gradus.isco(m) + 1e-2,
        #verbose = true,
        bins = bins,
        #maxrₑ = 500.0,
        maxrₑ = 30.0,
        callback = nothing,
        # resolution
        # numrₑ = 200,
        # #numrₑ = 50,
        # Nr = 3000,
        # abstol = 1e-10,
        # reltol = 1e-10,
        # callback = nothing,
        # kwargs...,
    )
    return f
end

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

function plot_all(data)
    incl_text = Printf.@sprintf " θ=%0.f h=%.1f" data.θ data.h
    p = plot(title = _format_metric(data.metric) * incl_text, legend = :topleft,  xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)}")
    line_colors = [:darkorange3, :mediumvioletred, :purple3, :blue3]
    for (i, (edd, f)) in enumerate(zip((0, 10, 20, 30), data.f))
        label_string = L"\dot{M} / \dot{M}_\textrm{edd} = " * string(edd / 100)
        plot!(p, data.bins, f, label = label_string, linecolor = line_colors[i])
    end
    p
end

bins = collect(range(0.1, 1.4, 200))

# models
model3 = LampPostModel(h = 3.0)
model6 = LampPostModel(h = 6.0)
model12 = LampPostModel(h = 12.0)

# inclinations
inclination15 = 15.0
inclination30 = 30.0
inclination60 = 60.0

#moving offset 0 -> 1 to avoid transfer functions' confusion
KWARGS = (; β₀ = 2)

################## i = 15 degrees

data_n1 = run_all_parameter_combinations(model3, inclination15, bins; KWARGS...)

pn1 = plot_all(data_n1)
display(pn1)

data_n2 = run_all_parameter_combinations(model6, inclination15, bins; KWARGS...)

pn2 = plot_all(data_n2)
display(pn2)

data_n3 = run_all_parameter_combinations(model12, inclination15, bins; KWARGS...)

pn3 = plot_all(data_n3)
display(pn3)

################## i = 30 degrees

data_k1 = run_all_parameter_combinations(model3, inclination30, bins; KWARGS...)

pk1 = plot_all(data_k1) 
display(pk1)


data_k2 = run_all_parameter_combinations(model6, inclination30, bins; KWARGS...)

pk2 = plot_all(data_k2) 
display(pk2)


data_k3 = run_all_parameter_combinations(model12, inclination30, bins; KWARGS...)

pk3 = plot_all(data_k3) 
display(pk3)

################## i = 60

data_p1 = run_all_parameter_combinations(model3, inclination60, bins; KWARGS...)

pp1 = plot_all(data_p1) 
display(pp1)

data_p2 = run_all_parameter_combinations(model6, inclination60, bins; KWARGS...)

pp2 = plot_all(data_p2) 
display(pp2)

data_p3 = run_all_parameter_combinations(model12, inclination60, bins; KWARGS...)

pp3 = plot_all(data_p3) 
display(pp3)

# put everything together for plotting
plot(pn3, pk3, pp3, pn2, pk2, pp2, pn1, pk1, pp1, layout = Plots.grid(3, 3), size = (1100, 1100))
