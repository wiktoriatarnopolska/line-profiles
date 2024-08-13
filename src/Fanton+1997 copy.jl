using Gradus
using Plots

#geometric thin disc
d = ThinDisc(0.0, 1e3)

m = KerrMetric(M=1.0, a = 0.998)


# define custom bins for g
bins = collect(range(0.1, 2.0, 300))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)

# inclination
x = SVector(0.0, 1000.0, deg2rad(30), 0.0)

model = LampPostModel(h = 10.0)

spectrum2 = PowerLawSpectrum(1.5)
spectrum3 = PowerLawSpectrum(2.0)
spectrum4 = PowerLawSpectrum(2.5)
spectrum5 = PowerLawSpectrum(3.0)



em_prof2 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum2,
    n_samples = 10_000,
    # sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

em_prof3 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum3,
    n_samples = 10_000,
    # sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

em_prof4 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum4,
    n_samples = 10_000,
    #sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

em_prof5 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum5,
    n_samples = 10_000,
    #sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

profile2 = Gradus.RadialDiscProfile(em_prof2)
profile3 = Gradus.RadialDiscProfile(em_prof3)
profile4 = Gradus.RadialDiscProfile(em_prof4)
profile5 = Gradus.RadialDiscProfile(em_prof5)


plot(profile2, label = "Γ = 1.5")
plot!(profile3, label = "Γ = 2.0")
plot!(profile4, label = "Γ = 2.5")
plot!(profile5, label = "Γ = 3.0")


em_prof = Gradus.emissivity_profile(
    m,
    d,
    model,
    #spectrum5,
    n_samples = 10_000,
    #sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

profile = Gradus.RadialDiscProfile(em_prof)


plot(profile3, label = "Γ = 2.0", 
title = "Test default vs. Γ = 2.0 vs. Γ = 3.0",
linewidth = 3)
plot!(profile, 
linestyle =:dash, label = "default (expectadly PowerLawSpectrum(2.0))")
plot!(profile5, label = "Γ = 3.0")