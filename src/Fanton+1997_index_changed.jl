# using Gradus
# using Plots
# using LaTeXStrings

# # inclination
# x = SVector(0.0, 1000.0, deg2rad(30), 0.0)


# #geometric thin disc
# d = ThinDisc(0.0, Inf)
# m = KerrMetric(M = 1.0, a = 0.0)
# # maximal integration radius
# maxrₑ = 20.0
# minrₑ = 6.0

# gs = range(0.0, 1.2, 500)


# # emissivity function
# ε1 =  ε(r) = r^(-2)
# _, flux1 = lineprofile(gs, ε1, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)

# ε2 = ε(r) = r^(-2.5)
# _, flux2 = lineprofile(gs, ε2, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)

# ε3 = ε(r) = r^(-3)
# _, flux3 = lineprofile(gs, ε3, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)


# # plot flux as a function of energy
# plot(gs, flux1, legend=true, label = "p = 2.0", colour =:hotpink2, title = "a = 0.0, r_in = 6.0, r_out = 20.0, θ = 30",  xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)}")
# plot!(gs, flux2, legend=true, colour =:tan1, label = "p = 2.5")
# plot!(gs, flux3, legend=true, colour =:darkorchid1, label = "p = 3.0")


# # inclination
# x = SVector(0.0, 1000.0, deg2rad(30), 0.0)


# #geometric thin disc
# d = ThinDisc(0.0, Inf)
# m = KerrMetric(M = 1.0, a = 0.998)
# # maximal integration radius
# maxrₑ = 20.0
# minrₑ = 1.23

# gs = range(0.0, 1.2, 500)


# # emissivity function
# ε1 =  ε(r) = r^(-2)
# _, flux1 = lineprofile(gs, ε1, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)

# ε2 = ε(r) = r^(-2.5)
# _, flux2 = lineprofile(gs, ε2, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)

# ε3 = ε(r) = r^(-3)
# _, flux3 = lineprofile(gs, ε3, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)


# # plot flux as a function of energy
# plot(gs, flux1, legend=true, label = "p = 2.0", colour =:hotpink2, title = "a = 0.998, r_in = 1.23, r_out = 20.0, θ = 30",  xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)}")
# plot!(gs, flux2, legend=true, colour =:tan1, label = "p = 2.5")
# plot!(gs, flux3, legend=true, colour =:darkorchid1, label = "p = 3.0")

using Gradus
using Plots
using LaTeXStrings


# define custom bins for g
bins = collect(range(0.1, 2.0, 300))

# define the plane to perform the binning over
#plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 20.0)

# inclination
x1 = SVector(0.0, 1000.0, deg2rad(5), 0.0)
x2 = SVector(0.0, 1000.0, deg2rad(30), 0.0)
x3 = SVector(0.0, 1000.0, deg2rad(45), 0.0)
x4 = SVector(0.0, 1000.0, deg2rad(85), 0.0)

#geometric thin disc
d = ThinDisc(0.0, Inf)

m = KerrMetric(M = 1.0, a = 0.998)
ε = ε(r) = r^(-2)
function calculate_line_profile(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        bins,
        ε,
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
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

function calculate_line_profile_binned(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        bins,
        ε,
        m,
        x,
        d;
        method = BinningMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        #plane=plane,
        minrₑ = 1.23,
        maxrₑ = 20.0,
        # resolution
        # numrₑ = 350,
        # Nr = 3000,
        # abstol = 1e-10,
        # reltol = 1e-10,
        # kwargs...,
    )
    return f
end


edrat1 = calculate_line_profile_binned(m, x1, d, bins)
edrat2 = calculate_line_profile_binned(m, x2, d, bins)
edrat3 = calculate_line_profile_binned(m, x3, d, bins)
edrat4 = calculate_line_profile_binned(m, x4, d, bins)

plot(bins, edrat1, label = "5 degrees", title = "Binning Method",  xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)}", colour =:darkorchid1,)
plot!(bins, edrat2, label = "30 degrees", colour =:tan1,)
plot!(bins, edrat3, label = "45 degrees", colour =:hotpink2,)
plot!(bins, edrat4, label = "85 degrees", colour =:midnightblue)

edrat1transfer = calculate_line_profile(m, x1, d, bins)
edrat2transfer = calculate_line_profile(m, x2, d, bins)
edrat3transfer = calculate_line_profile(m, x3, d, bins)
edrat4transfer = calculate_line_profile(m, x4, d, bins)

plot(bins, edrat1transfer, label = "5 degrees", title = "Transfer Functions Method",  xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)}", colour =:darkorchid1,)
plot!(bins, edrat2transfer, label = "30 degrees", colour =:tan1,)
plot!(bins, edrat3transfer, label = "45 degrees", colour =:hotpink2,)
plot!(bins, edrat4transfer, label = "85 degrees", colour =:midnightblue)


# #######################

# # inclination
# x = SVector(0.0, 1000.0, deg2rad(30), 0.0)

# #geometric thin disc
# d = ThinDisc(0.0, Inf)
# m = KerrMetric(M = 1.0, a = 0.998)

# minrₑ = 1.23

# gs = range(0.0, 1.2, 500)

# # emissivity function
# ε =  ε(r) = r^(-2)


# function calculate_line_profile(m, x, d, bins; kwargs...)
#     _, f = lineprofile(
#         bins,
#         ε,
#         m,
#         x,
#         d;
#         method = TransferFunctionMethod(),
#         #minrₑ = Gradus.isco(m) + 1e-2,
#         verbose = true,
#         minrₑ = minrₑ,
#         maxrₑ = 10.0,
#         # resolution
#         numrₑ = 350,
#         Nr = 3000,
#         abstol = 1e-10,
#         reltol = 1e-10,
#         kwargs...,
#     )
#     return f
# end


# edrat1 = calculate_line_profile(m, x, d, bins)


# function calculate_line_profile(m, x, d, bins; kwargs...)
#     _, f = lineprofile(
#         bins,
#         ε,
#         m,
#         x,
#         d;
#         method = TransferFunctionMethod(),
#         #minrₑ = Gradus.isco(m) + 1e-2,
#         verbose = true,
#         minrₑ = minrₑ,
#         maxrₑ = 20.0,
#         # resolution
#         numrₑ = 350,
#         Nr = 3000,
#         abstol = 1e-10,
#         reltol = 1e-10,
#         kwargs...,
#     )
#     return f
# end


# edrat2 = calculate_line_profile(m, x, d, bins)


# function calculate_line_profile(m, x, d, bins; kwargs...)
#     _, f = lineprofile(
#         bins,
#         ε,
#         m,
#         x,
#         d;
#         method = TransferFunctionMethod(),
#         #minrₑ = Gradus.isco(m) + 1e-2,
#         verbose = true,
#         minrₑ = minrₑ,
#         maxrₑ = 50.0,
#         # resolution
#         numrₑ = 350,
#         Nr = 3000,
#         abstol = 1e-10,
#         reltol = 1e-10,
#         kwargs...,
#     )
#     return f
# end


# edrat3 = calculate_line_profile(m, x, d, bins)


# function calculate_line_profile(m, x, d, bins; kwargs...)
#     _, f = lineprofile(
#         bins,
#         ε,
#         m,
#         x,
#         d;
#         method = TransferFunctionMethod(),
#         #minrₑ = Gradus.isco(m) + 1e-2,
#         verbose = true,
#         minrₑ = minrₑ,
#         maxrₑ = 100.0,
#         # resolution
#         numrₑ = 350,
#         Nr = 3000,
#         abstol = 1e-10,
#         reltol = 1e-10,
#         kwargs...,
#     )
#     return f
# end


# edrat4 = calculate_line_profile(m, x, d, bins)



# plot(bins, edrat1, label = "r_out = 10", title = "a = 0.998, r_in = 1.23, θ = 30",  xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\textrm{Flux~(arbitrary~units)}", colour =:darkorchid1,)
# plot!(bins, edrat2, label = "r_out = 20", colour =:tan1,)
# plot!(bins, edrat3, label = "r_out = 50", colour =:hotpink2,)
# plot!(bins, edrat4, label = "r_out = 100", colour =:midnightblue)
