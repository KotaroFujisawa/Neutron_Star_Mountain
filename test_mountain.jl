using Plots
using DifferentialEquations
using Dierckx
using QuadGK

include("Background_star.jl")
using .Background_star
include("Poisson_eq.jl")
using .Poisson_eq
include("Perturb_star.jl")
using .Perturb_star

function main()
    #simple NS model
    Γ = 2
    K = 46000     #cgs
    M0 = 1.989e33 #cgs
    G = 6.6743e-8 #cgs
    κ = 1.0e16    #cgs

    ℓ  = 2        #ℓ = m = 2 deformation
    β2 = ℓ*(ℓ+1)  #β^2 = ℓ(ℓ+1)
    β = sqrt(β)

    ρ0 = 2.0e15   #central density     r = r_min
    ρc = 2.0e14   #core-crust boundary r = r_c
    ρo = 1.0e11   #crust-ocian boundary r = r_o 
    ρ_min = 1.0e4 #surface desnity

    r_min = 0.1   #minimun radius
    r_c = 0.0     #core-crust boundary
    r_o = 0.0     #crust-ocian boundary
    R   = 0.0     #stellar radius
    r_max = 3.0e6 
    #make background star
    nr, R, M, r_g, ρ, p, cs2, m, dρ_dr, dp_dr, d2ρ_dr2 = (
        Background_star.make_bg_star_p(
            K, Γ, ρ0, ρ_min, r_min, r_max
        )
    )
    #r_g : radial grid

    #background star
    print("grid number = $nr, M = $(M/M0)  R = $(R/1.0e5) [km] \n")

    #interpolation by spline curve
    ρ_r = Spline1D(r_g, ρ)
 
    r_c = Background_star.cal_r_from_ρ(0.8*R,  ρc, ρ_r)
    r_o = Background_star.cal_r_from_ρ(0.95*R, ρo, ρ_r)

    print("r_c = $r_c  \n")
    print("r_o = $r_o  \n")

    #perturbation
    δρ = zeros(Float64, nr)
    δp = zeros(Float64, nr)
    δϕ = zeros(Float64, nr)
    dδϕ_dr = zeros(Float64, nr)

    #force A f_i = -Aρ ∇_i (r^2 Y_lm)
    A = 1e9
    fr(r) = -2A*r*ρ_r(r)
    ft(r) =  -A*r*ρ_r(r)

### fluid star ###

    δρ_r = Spline1D(r_g, δρ)

    max_ite = 50
    for ite=1:max_ite
        δρ = Perturb_star.δρ_fluid_star(r_g, ρ, δϕ, cs2, fr, ft, nr)
        δρ_r = Spline1D(r_g, δρ)
        δϕ, dδϕ_dr = Poisson_eq.Poisson(δρ_r, ℓ, r_min, R, r_g, nr)
        δρ_rc = δρ_r(r_c)
    end
    #ellipticity
    ε_f = Perturb_star.ellipticity(δρ_r, r_min, R)
    print("ε_f = $ε_f \n")

    plot(r_g, δρ)

#solid crust

    ξr = zeros(Float64, nr)
    ξt = zeros(Float64, nr)
    T1 = zeros(Float64, nr)
    T2 = zeros(Float64, nr)
    μ  = zeros(Float64, nr)
    dμ_dr = zeros(Float64, nr)

    for i=1:nr
#        if r_c < r_g[i] && r_g[i]< r_o
        μ[i] = κ*ρ[i]
#        end
    end

    μ_r     = Spline1D(r_g, μ)
    for i=1:nr 
        dμ_dr[i] = derivative(μ_r, r_g[i])
    end
    dμ_dr_r = Spline1D(r_g, dμ_dr)

    max_ite2 = 1
    for ite=1:max_ite2
        for i=1:nr
            if r_c < r_g[i] && r_g[i] < r_o
                δρ[i] = -(3ρ[i]/r_g[i] + dρ_dr[i])*ξr[i] + 3β*r_g[i]*ρ[i]/(2r_g[i])*ξt[i] - 3ρ[i]/(4μ[i])*T1[i] 
            else
                δρ[i] =  -(ρ[i]*δϕ[i] - ft(r_g[i])) / cs2[i]
            end
        end
    end
        
    plot(r_g, dμ_dr)

end

main()
