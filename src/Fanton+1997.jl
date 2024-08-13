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

m = KerrMetric(M = 1.0, a = 0.0)

model = LampPostModel(h = 10.0)

spectrum1 = PowerLawSpectrum(2.0)
spectrum2 = PowerLawSpectrum(2.5)
spectrum3 = PowerLawSpectrum(3.0)

em_prof1 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum1,
    n_samples = 10_000,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

em_prof2 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum2,
    n_samples = 10_000,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

em_prof3 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum3,
    n_samples = 10_000,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

prof1 = Gradus.RadialDiscProfile(em_prof1)
prof2 = Gradus.RadialDiscProfile(em_prof2)
prof3 = Gradus.RadialDiscProfile(em_prof3)

function calculate_line_profile(m, x, d, prof, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d,
        prof;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        minrₑ = 6.0,
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


edrat1 = calculate_line_profile(m, x, d, prof1, bins)
edrat2 = calculate_line_profile(m, x, d, prof2, bins)
edrat3 = calculate_line_profile(m, x, d, prof3, bins)


plot(bins, edrat1, label = "p = 2.0", title = "a = 0.0, r_in = 6.0, r_out = 20")
plot!(bins, edrat2, label = "p = 2.5")
plot!(bins, edrat3, label = "p = 3.0")


