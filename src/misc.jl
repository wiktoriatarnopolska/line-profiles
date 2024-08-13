using Gradus
using Plots
using Printf

_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.2f" m.a

_format_label(edd) = Printf.@sprintf "Ṁ / Ṁedd = %.1f" (edd / 100)

function em_prof(m, d, model; kwargs...)
    em_prof = Gradus.emissivity_profile(
        m,
        d,
        model,
        n_samples = 10_000,
        #sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
    )
    return em_prof
end


function run_all_parameter_combinations(m, model; kwargs...)

    @info "m = $m model = $(model)"


    # discs
    dthin = ThinDisc(Gradus.isco(m), Inf)
    d1 = ShakuraSunyaev(m, eddington_ratio = 0.1)
    d2 = ShakuraSunyaev(m, eddington_ratio = 0.2)
    d3 = ShakuraSunyaev(m, eddington_ratio = 0.3)

    @info "0%"
    em_prof0 = em_prof(m, dthin, model)
    prof0 = Gradus.RadialDiscProfile(em_prof0)

    @info "10%"
    em_prof10 = em_prof(m, d1, model)
    prof10 = Gradus.RadialDiscProfile(em_prof10)

    @info "20%"
    em_prof20 = em_prof(m, d2, model)
    prof20 = Gradus.RadialDiscProfile(em_prof20)

    @info "30%"
    em_prof30 = em_prof(m, d3, model)
    prof30 = Gradus.RadialDiscProfile(em_prof30)

    return (; metric = m, f = [prof0, prof10, prof20, prof30], h = model.h)

end

function plot_all(data)
    incl_text = Printf.@sprintf "h=%0.f" data.h
    p = plot(title = _format_metric(data.metric) * incl_text, legend = :bottomleft, xlims = (1.0, 30.0))
    for (edd, f) in zip((0, 10, 20, 30), data.f)
        plot!(p, f, label = _format_label(edd))
    end
    p
end



# models
model3 = LampPostModel(h = 3.0)
model6 = LampPostModel(h = 6.0)
model12 = LampPostModel(h = 12.0)

################## a = 0.0


m_n1 = KerrMetric(M = 1.0, a = 0.0)
data_n1 = run_all_parameter_combinations(m_n1, model3)

pn1 = plot_all(data_n1)
display(pn1)

m_n2 = KerrMetric(M = 1.0, a = 0.0)
data_n2 = run_all_parameter_combinations(m_n2, model6)


pn2 = plot_all(data_n2)
display(pn2)

m_n3 = KerrMetric(M = 1.0, a = 0.0)
data_n3 = run_all_parameter_combinations(m_n3, model12)

pn3 = plot_all(data_n3)
display(pn3)

################## a = 0.9

m_k1 = KerrMetric(M = 1.0, a = 0.9)

data_k1 = run_all_parameter_combinations(m_k1, model3)

pk1 = plot_all(data_k1) 
display(pk1)

m_k2 = KerrMetric(M = 1.0, a = 0.9)
data_k2 = run_all_parameter_combinations(m_k2, model6)

pk2 = plot_all(data_k2) 
display(pk2)

m_k3 = KerrMetric(M = 1.0, a = 0.9)
data_k3 = run_all_parameter_combinations(m_k3, model12)

pk3 = plot_all(data_k3) 
display(pk3)

################## a = 0.99

m_p1 = KerrMetric(M = 1.0, a = 0.99)
data_p1 = run_all_parameter_combinations(m_p1, model3)

pp1 = plot_all(data_p1) 
display(pp1)

m_p2 = KerrMetric(M = 1.0, a = 0.99)
data_p2 = run_all_parameter_combinations(m_p2, model6)

pp2 = plot_all(data_p2) 
display(pp2)

m_p3 = KerrMetric(M = 1.0, a = 0.99)
data_p3 = run_all_parameter_combinations(m_p3, model12)

pp3 = plot_all(data_p3) 
display(pp3)

# put everything together
plot(pn3, pk3, pp3, pn2, pk2, pp2, pn1, pk1, pp1,layout = grid(3, 3), size = (1100, 1100))

