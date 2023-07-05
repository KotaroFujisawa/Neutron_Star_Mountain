using Plots
using DifferentialEquations
using Dierckx
using QuadGK
using PyCall
@pyimport matplotlib.pyplot as plt

include("Background_star.jl")
using .Background_star
include("Poisson_eq.jl")
using .Poisson_eq
include("Perturb_star.jl")
using .Perturb_star

function main()
    #simple NS model
    Γ = 2
    K = 42500     #cgs
    M0 = 1.989e33 #cgs solar mass
    G = 6.6743e-8 #cgs
    κ = 1.0e16    #cgs

    ℓ  = 2        #ℓ = m = 2 deformation
    β2 = ℓ*(ℓ+1)  #β^2 = ℓ(ℓ+1)
    β = sqrt(β2)

    ρ0 = 2.185e15   #central density     r = r_min
    ρc = 2.0e14   #core-crust boundary r = r_c
    ρo = 2.18e8 #2.18e8   #crust-ocian boundary r = r_o 
    ρ_min = 1.0e4 #surface desnity

    r_min = 0.1   #minimun radius
    r_c = 0.0     #core-crust boundary
    r_o = 0.0     #crust-ocian boundary
    R   = 0.0     #stellar radius
    r_max = 3.0e6 

    #make background star (N=1 polytropic star)
    nr, R, M, r_g, ρ, p, cs2, m, dρ_dr, dp_dr, d2ρ_dr2 = (
        Background_star.make_bg_star_p(
            K, Γ, ρ0, ρ_min, r_min, r_max
        )
    )
    #r_g : radial grid

    #physical quantities of background star
    print("grid number = $nr, M = $(M/M0)  R = $(R/1.0e5) [km] \n")

    #interpolation by spline curve
    ρ_r     = Spline1D(r_g, ρ)
    dρ_dr_r = Spline1D(r_g, dρ_dr)    #dρ / dr
    d2ρ_dr2_r= Spline1D(r_g, d2ρ_dr2) #d^2ρ / dr^2
    cs2_r   = Spline1D(r_g, cs2)  

    dcs2_dr = zeros(Float64, nr)
    for i=1:nr 
        dcs2_dr[i] = derivative(cs2_r, r_g[i])
    end
    dcs2_dr_r = Spline1D(r_g, dcs2_dr)

    #calculate r_o & r_c
    r_c = Background_star.cal_r_from_ρ(0.8*R,  ρc, ρ_r)
    r_o = Background_star.cal_r_from_ρ(0.95*R, ρo, ρ_r)

    print("r_c = $r_c  \n")
    print("r_o = $r_o  \n")

    out = open("rho.dat","w")
    for i=1:nr
        print(out, "$(r_g[i])  $(ρ[i]) \n")
    end
    close(out)

    
    #perturbation
    δρ = zeros(Float64, nr)
    δp = zeros(Float64, nr)
    δϕ = zeros(Float64, nr)
    dδϕ_dr = zeros(Float64, nr)

    #force A f_i = -Aρ ∇_i (r^2 Y_lm)
    A = 0.74e7
    fr(r) =-2A*r*ρ_r(r)
    ft(r) =-A*r*ρ_r(r)
    
    #force B f_i = Bρr ∇_i Y_lm 
    # fr = 0.0
    # ft = Bρ 
#    B = 2*10.2e9/1.5068 /sqrt(0.979)
#    fr(r) = 0.0
#    ft(r) = B*ρ_r(r)

### fluid star ###

    δρ_r = Spline1D(r_g, δρ)

    max_ite = 50
    for ite=1:max_ite
        δρ = Perturb_star.δρ_fluid_star(r_g, ρ, δϕ, cs2, fr, ft, nr)
        δρ_r = Spline1D(r_g, δρ)
        δϕ, dδϕ_dr = Poisson_eq.Poisson_BVP(δρ_r, ℓ, r_min, R, r_g, nr)
        δρ_rc = δρ_r(r_c)
    end
    #ellipticity
    ε_f = Perturb_star.ellipticity(δρ_r, r_min, R)


#solid crust

    ξr = zeros(Float64, nr)
    ξt = zeros(Float64, nr)
    T1 = zeros(Float64, nr)
    T2 = zeros(Float64, nr)
    μ  = zeros(Float64, nr)
    dμ_dr = zeros(Float64, nr)
 
    for i=1:nr
        μ[i] = κ*ρ[i]
    end

    μ_r     = Spline1D(r_g, μ)
    for i=1:nr 
        dμ_dr[i]   = derivative(μ_r, r_g[i])
    end
    dμ_dr_r   = Spline1D(r_g, dμ_dr)

    max_ite2 = 15

    for ite=1:max_ite2
        for i=1:nr
            if r_c < r_g[i] && r_g[i] < r_o
                δρ[i] = -(3ρ[i]/r_g[i] + dρ_dr[i])*ξr[i] + 3β*ρ[i]/(2r_g[i])*ξt[i] - 3ρ[i]/(4μ[i])*T1[i] 
            else
                δρ[i] =  -(ρ[i]*δϕ[i] - ft(r_g[i])*r_g[i]) / cs2[i]
            end
        end
        δρ_r     = Spline1D(r_g, δρ)
        δϕ_r     = Spline1D(r_g, δϕ)
        dδϕ_dr_r = Spline1D(r_g, dδϕ_dr)

        T1, T2, ξr, ξt = Perturb_star.cal_ξ_T_BVP(ρ_r, dρ_dr_r, d2ρ_dr2_r, cs2_r, dcs2_dr_r, 
        μ_r, dμ_dr_r, δϕ_r, dδϕ_dr_r, fr, ft, β2, r_c, r_o, r_g, nr)

        δϕ, dδϕ_dr = Poisson_eq.Poisson_BVP(δρ_r, ℓ, r_min, R, r_g, nr)

    end

    #δp
    for i=1:nr
        δp[i] = cs2[i] * δρ[i]
    end

    ε_s = Perturb_star.ellipticity(δρ_r, r_min, R)
    println("ε_f = $ε_f")
    println("ε_s = $ε_s")
    println("|ε_s - ε_f| = ", abs(ε_s - ε_f))


    nθ = 100
    nφ = 100


    θ = range(0.0, stop=π,  length=nθ)
    φ = range(0.0, stop=2π, length=nφ) 

    #strain
    σ2 = zeros(Float64, nr, nθ, nφ)
    σ  = zeros(Float64, nr)
    σr1 = zeros(Float64, nr)
    σr2 = zeros(Float64, nr)
    σr3 = zeros(Float64, nr)
    for i=1:nr
        for j=1:nθ
            for k = 1:nφ
                if r_c ≤ r_g[i] && r_g[i] ≤ r_o
                    σ2[i,j,k] = 5/(256π) * (6*sin(θ[j])^2*
                        (3*sin(θ[j])^2*cos(2φ[k])^2 * (T1[i]/μ[i])^2
                        +4*(3+cos(2θ[j]) - 2*sin(θ[j])^2*cos(4φ[k]))*(T2[i]/μ[i])^2
                        )  + (35+28*cos(2θ[j]) + cos(4θ[j]) 
                        + 8*sin(θ[j])^4*cos(4φ[k]))*(ξt[i]/r_g[i])^2 
                    )
                else
                    σ2[i,j,k] = 0.0
                end
            end 
        end
    end

    for i=1:nr
        if r_c ≤ r_g[i] && r_g[i] ≤ r_o
            σr1[i] = sqrt(45.0/128π) * abs(T1[i] / μ[i])
            σr2[i] = sqrt(1215/512π) * abs(T2[i] / μ[i])
            σr3[i] = sqrt(5.0/4π) * abs(ξt[i]/ r_g[i])
        else
            σr1[i] = 0.0
            σr2[i] = 0.0
            σr3[i] = 0.0
        end
    end

    out = open("sigma123.dat","w")
    for i=1:nr
        print(out, "$(r_g[i])  $(σr1[i])  $(σr2[i]) $(σr3[i]) \n")
    end
    close(out)

    out = open("T1T2.dat","w")
    for i=1:nr
        tmp_T1 =( ((ρ[i]*δϕ[i] - ft(r_g[i])*r_g[i]) - 
        cs2[i]*( (3ρ[i]/r_g[i] + dρ_dr[i])*ξr[i] - 3β*ρ[i]/2r_g[i]*ξt[i]))
        / (1 + 3cs2[i]*ρ[i]/(4μ[i])))
        print(out, "$(r_g[i])  $(δp[i] - T1[i]) $(δp[i])  $(T1[i])  $(tmp_T1)  $(T2[i]) \n")
    end
    close(out)

    σ2_max = maximum(σ2)
    #maximum strain
    println("σ_max = ", sqrt(σ2_max))

    for i=1:nr
        σ[i] = sqrt(maximum(σ2[i,:,:]))
    end

    #plot by using matplotlib

    fig, ax = plt.subplots()
    ax.plot(r_g, T1)
    ax.grid(true)
    ax.set_xlabel("r [cm]")
    ax.set_xlim(r_c, r_o)
    ax.set_ylabel("T1")
    plt.savefig("T1.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, δp - T1)
    ax.grid(true)
    ax.set_ylabel("δp - T1")
    ax.set_xlim(0.0, R)
    plt.savefig("deltap_T1.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, T2)
    ax.grid(true)
    ax.set_ylabel("T2")
    ax.set_xlim(r_c, r_o)
    plt.savefig("T2.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, ξr)
    ax.grid(true)
    ax.set_ylabel("ξr")
    ax.set_xlim(r_c, r_o)
    plt.savefig("xi_r.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, ξt)
    ax.grid(true)
    ax.set_ylabel("ξt")
    ax.set_xlim(r_c, r_o)
    plt.savefig("xi_t.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, δρ)
    ax.grid(true)
    ax.set_ylabel("δρ")
    ax.set_xlim(r_min, R)
    plt.savefig("delta_rho.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, δϕ)
    ax.grid(true)
    ax.set_ylabel("δϕ")
    ax.set_xlim(r_min, R)
    plt.savefig("delta_phi.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, σ)
    ax.grid(true)
    ax.set_ylabel("|σ|")
    ax.set_xlim(0.99r_o, r_o)
    plt.savefig("sigma.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.plot(r_g, σr1)
    ax.plot(r_g, σr2)
    ax.plot(r_g, σr3)
    ax.grid(true)
    ax.set_ylabel("|σ|")
    ax.set_xlim(0.99r_o, r_o)
    ax.set_ylim(1.0e-6, 0.1)
    plt.savefig("sigma123.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.plot(r_g, ρ)
    ax.grid(true)
    ax.set_ylabel("ρ")
    ax.set_xlim(0.999r_o, 1.001r_o)
    plt.savefig("density.pdf")
    plt.close()

end

main()