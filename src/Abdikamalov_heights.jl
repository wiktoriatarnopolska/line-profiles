using Gradus
using Plots
using Printf
using LaTeXStrings
using Plots.PlotMeasures

_format_metric(m::JohannsenMetric) = Printf.@sprintf "α13=%.2f, ϵ3=%.1f" m.α13 m.ϵ3
_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.2f" m.a

function calculate_line_profile(m, x, d, prof, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d,
        prof;
        method = TransferFunctionMethod(),
        minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        maxrₑ = 500.0,
        # resolution
        numrₑ = 200,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
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



function run_all_parameter_combinations(m, θ, model, bins; kwargs...)

    x = SVector(0.0, 1000.0, deg2rad(θ), 0.0)

    spectrum = PowerLawSpectrum(2.0)

    @info "m = $m θ = $(θ)"

    # discs
    dthin = ThinDisc(0.0, Inf)
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
    #incl_text = Printf.@sprintf "h=%.1f" data.h
    p = plot(title =  _format_metric(data.metric), legend = :topleft, xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)} ")
    line_colors = [:darkorange3, :mediumvioletred, :purple3, :blue3]
    for (i, (edd, f)) in enumerate(zip((0, 10, 20, 30), data.f))
        label_string = L"\dot{M} / \dot{M}_\textrm{edd} = " * string(edd / 100)
        plot!(p, data.bins, f, label = label_string, linecolor = line_colors[i])
    end
    p
end

# define custom bins for g
bins = collect(range(0.1, 1.4, 200))

INCLINATION = 70.0
# check for β = 2
KWARGS = (; β₀ = 2)

m_n1 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.5, ϵ3 = 1.0)
model5 = LampPostModel(h = 6.0)
data_n1 = run_all_parameter_combinations(m_n1, INCLINATION, model5, bins; KWARGS...)

pn1 = plot_all(data_n1)
display(pn1)

m_n2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.5, ϵ3 = 0.5)
model10 = LampPostModel(h = 6.0)
data_n2 = run_all_parameter_combinations(m_n2, INCLINATION, model10, bins; KWARGS...)


pn2 = plot_all(data_n2)
display(pn2)

m_n3 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.5, ϵ3 = -1.0)
model15 = LampPostModel(h = 6.0)
data_n3 = run_all_parameter_combinations(m_n3, INCLINATION, model15, bins; KWARGS...)

pn3 = plot_all(data_n3)
display(pn3)


# put everything together
#plot(pn1, pn2, pn3, layout = grid(1, 3, heights=[0.3]), size = (1200, 1200), left_margin = [5mm 0mm], right_margin = [5mm 0mm])
plot(pn1, pn2, pn3, layout = grid(3, 1), size = (600, 1200), left_margin = [12mm 0mm])

m_p1 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -1.0, ϵ3 = 0.0)
model5 = LampPostModel(h = 6.0)
data_p1 = run_all_parameter_combinations(m_p1, INCLINATION, model5, bins; KWARGS...)

pp1 = plot_all(data_p1)
display(pp1)

m_p2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.0, ϵ3 = 0.0)
model10 = LampPostModel(h = 6.0)
data_p2 = run_all_parameter_combinations(m_p2, INCLINATION, model10, bins; KWARGS...)


pp2 = plot_all(data_p2)
display(pp2)

m_p3 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 1.0, ϵ3 = 0.0)
model15 = LampPostModel(h = 6.0)
data_p3 = run_all_parameter_combinations(m_p3, INCLINATION, model15, bins; KWARGS...)

pp3 = plot_all(data_p3)
display(pp3)


# put everything together
#plot(pn1, pn2, pn3, layout = grid(1, 3, heights=[0.3]), size = (1200, 1200), left_margin = [5mm 0mm], right_margin = [5mm 0mm])
plot(pp1, pp2, pp3, layout = grid(3, 1), size = (600, 1200), left_margin = [12mm 0mm])

m_k1 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -1.0, ϵ3 = 1.0)
model5 = LampPostModel(h = 6.0)
data_k1 = run_all_parameter_combinations(m_k1, INCLINATION, model5, bins; KWARGS...)

pk1 = plot_all(data_k1)
display(pk1)

m_k2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.5, ϵ3 = -1.0)
model10 = LampPostModel(h = 6.0)
data_k2 = run_all_parameter_combinations(m_k2, INCLINATION, model10, bins; KWARGS...)


pk2 = plot_all(data_k2)
display(pk2)

m_k3 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 1.0, ϵ3 = -1.0)
model15 = LampPostModel(h = 6.0)
data_k3 = run_all_parameter_combinations(m_k3, INCLINATION, model15, bins; KWARGS...)

pk3 = plot_all(data_k3)
display(pk3)


# put everything together
#plot(pn1, pn2, pn3, layout = grid(1, 3, heights=[0.3]), size = (1200, 1200), left_margin = [5mm 0mm], right_margin = [5mm 0mm])
plot(pk1, pk2, pk3, layout = grid(3, 1), size = (600, 1200), left_margin = [12mm 0mm])

plot(pn1, pk1, pp3, pn2, pp2, pk2, pn3, pk3, pp1, layout = grid(3, 3), size = (1100, 1100), left_margin = [5mm 0mm])

#========================================================#

using Gradus
using Plots
using Printf
using LaTeXStrings
using Plots.PlotMeasures

_format_metric(m::JohannsenMetric) = Printf.@sprintf "a=%.2f, α13=%.2f" m.a m.α13
_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.2f" m.a

function calculate_line_profile(m, x, d, prof, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d,
        prof;
        method = TransferFunctionMethod(),
        minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        maxrₑ = 500.0,
        # resolution
        numrₑ = 200,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
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



function run_all_parameter_combinations(m, θ, model, bins; kwargs...)

    x = SVector(0.0, 1000.0, deg2rad(θ), 0.0)

    spectrum = PowerLawSpectrum(2.0)

    @info "m = $m θ = $(θ)"

    # discs
    dthin = ThinDisc(0.0, Inf)
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
    incl_text = Printf.@sprintf "h=%.1f" data.h
    p = plot(title = incl_text, legend = :topleft, xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)} ")
    line_colors = [:darkorange3, :mediumvioletred, :purple3, :blue3]
    for (i, (edd, f)) in enumerate(zip((0, 10, 20, 30), data.f))
        label_string = L"\dot{M} / \dot{M}_\textrm{edd} = " * string(edd / 100)
        plot!(p, data.bins, f, label = label_string, linecolor = line_colors[i])
    end
    p
end

# define custom bins for g
bins = collect(range(0.1, 1.4, 200))

INCLINATION = 70.0
# check for β = 2
KWARGS = (; β₀ = 2)

m_n1 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
model5 = LampPostModel(h = 5.0)
data_n1 = run_all_parameter_combinations(m_n1, INCLINATION, model5, bins; KWARGS...)

pn1 = plot_all(data_n1)
display(pn1)

m_n2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
model10 = LampPostModel(h = 10.0)
data_n2 = run_all_parameter_combinations(m_n2, INCLINATION, model10, bins; KWARGS...)


pn2 = plot_all(data_n2)
display(pn2)

m_n3 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
model15 = LampPostModel(h = 15.0)
data_n3 = run_all_parameter_combinations(m_n3, INCLINATION, model15, bins; KWARGS...)

pn3 = plot_all(data_n3)
display(pn3)


# put everything together
#plot(pn1, pn2, pn3, layout = grid(1, 3, heights=[0.3]), size = (1200, 1200), left_margin = [5mm 0mm], right_margin = [5mm 0mm])
plot(pn1, pn2, pn3, layout = grid(3, 1), size = (600, 1200), left_margin = [12mm 0mm])
