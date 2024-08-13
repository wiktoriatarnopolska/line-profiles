using Gradus
using Plots
using Test



m = KerrMetric(1.0, 0.998)
d = ThinDisc(0.0, 100.0)

model = LampPostModel(h = 10.0)

# #calculating emissivity profile without defining the spectrum to see if it uses the default Γ value (Γ = 2)
# em_prof1 = Gradus.emissivity_profile(
# m, 
# d, 
# model, 
# n_samples = 100_000, 
# sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
# )

# #first profile, default photon index Γ
# profile1 = Gradus.RadialDiscProfile(em_prof1)

# #defining the spectrum
# spectrum = PowerLawSpectrum(2.0)

# #secondly, passing spectrum as an argument 
# em_prof2 = Gradus.emissivity_profile(
# m, 
# d, 
# model, 
# spectrum,
# n_samples = 100_000, 
# sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
# )

# plot(profile1, label = "default photon iondex", linecolor = "orchid3", legend = true, title = "Test for the photon index")
# plot!(profile2, label = "passing spec argument", linecolor = "royalblue4", legend = true)

# #profile number 2 with defined Γ = 2
# profile2 = Gradus.RadialDiscProfile(em_prof2)

# #testing if they are equal
# @test profile1.ε ≈ profile2.ε


#======================================================#

spectrum1 = PowerLawSpectrum(1.5)

em_prof1 = Gradus.emissivity_profile(
m, 
d, 
model, 
spectrum1,
n_samples = 500_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

profile1 = Gradus.RadialDiscProfile(em_prof1)

spectrum2 = PowerLawSpectrum(2.0)

em_prof2 = Gradus.emissivity_profile(
m, 
d, 
model, 
spectrum2,
n_samples = 500_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

profile2 = Gradus.RadialDiscProfile(em_prof2)

spectrum3 = PowerLawSpectrum(2.5)

em_prof3 = Gradus.emissivity_profile(
m, 
d, 
model, 
spectrum3,
n_samples = 500_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

profile3 = Gradus.RadialDiscProfile(em_prof3)

spectrum4 = PowerLawSpectrum(3.0)

em_prof4 = Gradus.emissivity_profile(
m, 
d, 
model, 
spectrum4,
n_samples = 500_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

profile4 = Gradus.RadialDiscProfile(em_prof4)


plot(profile1, label = "Γ = 1.5", linecolor = "darkorange1", legend = true, xlabel = "Radius (rg)", ylabel = "ϵ(r) (arbitrary units)")
plot!(profile2, label = "Γ = 2.0", linecolor = "slateblue1", legend = true)
plot!(profile3, label = "Γ = 2.5", linecolor = "magenta3", legend = true)
plot!(profile4, label = "Γ = 3.0", linecolor = "purple4", legend = true)

m = KerrMetric(1.0, 0.998)
d = ThinDisc(0.0, 100.0)

model = LampPostModel(h = 10.0)
#calculating emissivity profile without defining the spectrum to see if it uses the default Γ value (Γ = 2)
em_prof1 = Gradus.emissivity_profile(
m, 
d, 
model, 
n_samples = 500_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

#first profile, default photon index Γ
profile1 = Gradus.RadialDiscProfile(em_prof1)

#defining the spectrum
spectrum = PowerLawSpectrum(2.0)

#secondly, passing spectrum as an argument 
em_prof2 = Gradus.emissivity_profile(
m, 
d, 
model, 
spectrum,
n_samples = 500_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

#profile number 2 with defined Γ = 2
profile2 = Gradus.RadialDiscProfile(em_prof2)

plot(profile1, label = "default photon iondex", linecolor = "plum", legend = true, xlabel = "Radius (rg)", ylabel = "ϵ(r) (arbitrary units)", linewidth = 5)
plot!(profile2, label = "passing spec argument", linecolor = "purple", legend = true)


#testing if they are equal
@test profile1.ε ≈ profile2.ε