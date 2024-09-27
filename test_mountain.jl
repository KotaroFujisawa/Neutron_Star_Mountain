using BenchmarkTools
using Plots
using DifferentialEquations
using Dierckx
using QuadGK
using Printf
using PyCall
#using Conda
#Conda.add("matplotlib")
@pyimport matplotlib.pyplot as plt

include("Background_star.jl")
using .Background_star
include("Poisson_eq.jl")
using .Poisson_eq
include("Perturb_star.jl")
using .Perturb_star
include("Force_density.jl")
using .Force_density
function main()
    # N = 1 polytropic NS model
    Γ = 2
    M0 = 1.989e33 #cgs solar mass
    G = 6.6743e-8 #cgs
    κ = 1.0e16    #cgs
    radius = 1.000000e6 #cgs

    ρ_ave = 1.4*M0 / (4*π/3*radius^3) # average density
    ρ0 = ρ_ave  / (3/π / π )
    K = 4π*G / (2) * radius^2 / π^2
#
#    print(ρ0, " ", K)
#    return 0
#    ρ0 = 2.185e15 #central density at r = r_min
#    K = 42500     #cgs

    ℓ  = 2        #ℓ = m = 2 deformation
    β2 = ℓ*(ℓ+1)  #β^2 = ℓ(ℓ+1)
    β = sqrt(β2)

    ρc = 2.0e14   #core-crust boundary at r = r_c
    ρo = 2.18e8   #crust-ocian boundary at r = r_o 
    ρ_min = 1.0e4 #surface desnity

    r_min = 0.1   #minimun radius
    r_c = 0.0     #core-crust boundary
    r_o = 0.0     #crust-ocian boundary
    R   = 0.0     #stellar radius
    r_max = 3.0e6 #maximum radius of the domain

    # mesh size (accuracy) for background star default:1.0e-12 
    reltol = 1.0e-12
    abstol = 1.0e-12

    # convergence criterion for main loop
    eps_c = 1.0e-10

    ε_f = 0.0     # ellpiticity for fluid star
    ε_s = 0.0     # ellpiticity for solid star
    # type of force density
    type_force = 1
    # type of bc 1:normal 3: surface current model
    type_BC = 1
    # coefficients for surface current j_0: core-crsut, j_1: crust-ocean
    # coefficients for solenoidal-irrotational force j_0: irr, j_1: sol 
    j_0 =  0.0e0
    j_1 =  1.0e0
    # make background star (N=1 polytropic star)
    nr, R, M, r_g, ρ, p, cs2, m, dρ_dr, dp_dr, d2ρ_dr2 = (
        Background_star.make_bg_star_p(
            K, Γ, ρ0, ρ_min, r_min, r_max,  reltol, abstol
        )
    )

    # physical quantities of background star
    print("grid number = $nr, M = $(M/M0)  R = $(R/1.0e5) [km] \n")

    # interpolation using spline curve
    ρ_r     = Spline1D(r_g, ρ)
    dρ_dr_r = Spline1D(r_g, dρ_dr)    #dρ / dr
    d2ρ_dr2_r= Spline1D(r_g, d2ρ_dr2) #d^2ρ / dr^2
    cs2_r   = Spline1D(r_g, cs2)  

    dcs2_dr = zeros(Float64, nr)
    for i=1:nr 
        dcs2_dr[i] = derivative(cs2_r, r_g[i])
    end
    dcs2_dr_r = Spline1D(r_g, dcs2_dr)

    # calculate r_o at ρc & r_c at ρo
    r_c = Background_star.cal_r_from_ρ(0.8*R,  ρc, ρ_r)
    r_o = Background_star.cal_r_from_ρ(0.95*R, ρo, ρ_r)

    print("r_c = $r_c,  r_o = $r_o \n")

    # quantities for perturbation
    δρ = zeros(Float64, nr)
    δρ_f = zeros(Float64, nr)
    δp = zeros(Float64, nr)
    δϕ_f = zeros(Float64, nr)
    δϕ = zeros(Float64, nr)
    dδϕ_dr = zeros(Float64, nr)

    ξr = zeros(Float64, nr)
    ξt = zeros(Float64, nr)
    T1 = zeros(Float64, nr)
    T2 = zeros(Float64, nr)
    μ  = zeros(Float64, nr)
    dμ_dr = zeros(Float64, nr)

    fr_g = zeros(Float64, nr)
    ft_g = zeros(Float64, nr)
    # radial functions of magnetic stress tensor 
    M1 = zeros(Float64, nr)
    M2 = zeros(Float64, nr)

    σ  = zeros(Float64, nr)
    σr1 = zeros(Float64, nr)
    σr2 = zeros(Float64, nr)
    σr3 = zeros(Float64, nr)


    for i=1:nr
        μ[i] = κ*ρ[i]
    end
    μ_r     = Spline1D(r_g, μ)
    for i=1:nr 
        dμ_dr[i]   = derivative(μ_r, r_g[i])
    end
    dμ_dr_r   = Spline1D(r_g, dμ_dr)

    #coefficient to normalize the force
    A = 0.0

    # set force 
    A, fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc = (
        Force_density.set_force_density(type_force, ρ_r, dρ_dr_r, r_c, r_o, R, A, j_0, j_1)
    )
    m_rr_co, m_rr_cr, m_rr_oc, m_rth_co, m_rth_cr, m_rth_oc = (
        Force_density.set_magnetic_stress(A, r_c, r_o, R, j_0, j_1)
    )

    out = open("force.dat","w")
    for i=1:nr
        if r_g[i] < r_c  # core
            print(out, "$(r_g[i])  $(fr_co(r_g[i]))  $(ft_co(r_g[i]))\n")
        elseif r_g[i] <= r_o # crust
            print(out, "$(r_g[i])  $(fr_cr(r_g[i]))  $(ft_cr(r_g[i]))\n")
        else # ocean 
            print(out, "$(r_g[i])  $(fr_oc(r_g[i]))  $(ft_oc(r_g[i]))\n")
        end
    end
    close(out)

    # iteration
    global_ite = 0
    max_global_ite = 30
    max_ite1 = 100
    max_ite2 = 100

    for global_ite=1:max_global_ite
    ### fluid star ###

        dδϕ_dr_o = 0.0
        dδϕ_dr_n = 0.0

        δρ_r = Spline1D(r_g, δρ_f)
        guess = 0.0
        for ite=1:max_ite1
            δρ_f = Perturb_star.δρ_fluid_star(r_c, r_o, r_g, ρ, δϕ_f, cs2, 
                fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc, nr)
            δρ_r = Spline1D(r_g, δρ_f)
            δϕ_f, dδϕ_dr = Poisson_eq.Poisson_BVP0(δρ_r, ℓ, r_min, R, r_g, nr, guess)
            dδϕ_dr_n = dδϕ_dr[1]
            diff_dδϕ_dr = abs(dδϕ_dr_n - dδϕ_dr_o) / abs(dδϕ_dr_n)
            if(ite > 3)
                guess = (0.5*dδϕ_dr_n + 0.5*dδϕ_dr_o) 
            end
            println("diff = $diff_dδϕ_dr", " $dδϕ_dr_n $dδϕ_dr_o", "  guess = $guess")
            if(diff_dδϕ_dr < 1.0e-10) 
                println("Converge!")
                break
            end
            dδϕ_dr_o = dδϕ_dr_n
            ε_f = Perturb_star.ellipticity(δρ_r, r_min, R)
            println("ε_f = $ε_f")
    
        end
        #
        ε_f = Perturb_star.ellipticity(δρ_r, r_min, R)
    #solid crust

        guess = 0.0
        guess2 = [0.0, 0.0, 0.0]
        init2_o = [0.0, 0.0, 0.0]
        δρ_o = zeros(Float64, nr)
        δϕ, dδϕ_dr = Poisson_eq.Poisson_BVP0(δρ_r, ℓ, r_min, R, r_g, nr, guess)
        for ite=1:max_ite2
            println("ite = ", ite)

            δρ = Perturb_star.δρ_solid_star(r_c, r_o, r_g, ρ, dρ_dr, δϕ, cs2, 
                fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc, ξr, ξt, T1, μ, β, nr)

            δρ_r     = Spline1D(r_g, δρ)
            δϕ_r     = Spline1D(r_g, δϕ)
            dδϕ_dr_r = Spline1D(r_g, dδϕ_dr)

            T1, T2, ξr, ξt, init2 = Perturb_star.cal_ξ_T_BVP0!(ρ_r, dρ_dr_r, d2ρ_dr2_r, cs2_r, dcs2_dr_r, 
                μ_r, dμ_dr_r, δϕ_r, dδϕ_dr_r, 
                fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc, 
                m_rr_co, m_rr_cr, m_rr_oc, m_rth_co, m_rth_cr, m_rth_oc,
                β2, r_c, r_o, r_g, nr, guess2, type_BC)

            δϕ, dδϕ_dr = Poisson_eq.Poisson_BVP0(δρ_r, ℓ, r_min, R, r_g, nr, guess)
            dδϕ_dr_n = dδϕ_dr[1]
            diff_dδϕ_dr = abs(dδϕ_dr_n - dδϕ_dr_o) / abs(dδϕ_dr_n)
            if(ite > 1)
                guess = (0.5*dδϕ_dr_n + 0.5*dδϕ_dr_o) 
                guess2 = (0.5*init2 + 0.5*init2_o)
            end
            println("diff = $diff_dδϕ_dr", " $dδϕ_dr_n $dδϕ_dr_o", "  guess = $guess")
            dδϕ_dr_o =  dδϕ_dr_n
            init2_o = init2
            if(diff_dδϕ_dr < eps_c) 
                println("Converge!")
                break
            end
            dδϕ_dr_o = dδϕ_dr_n 
            δρ_o = δρ
            ε_s = Perturb_star.ellipticity(δρ_r, r_min, R)
            println("ε_f = $ε_f", " ε_s = $ε_s", " |ε_s - ε_f| = ", abs(ε_s - ε_f))
        end

        #δp
        for i=1:nr
            δp[i] = cs2[i] * δρ[i]
        end

        ε_s = Perturb_star.ellipticity(δρ_r, r_min, R)
        println("ε_f = $ε_f")
        println("ε_s = $ε_s")
        println("|ε_s - ε_f| = ", abs(ε_s - ε_f))

        σ2_max = Perturb_star.cal_σ!(σ, σr1, σr2, σr3, r_c, r_o, r_g, T1, T2, μ, ξr, ξt, nr)

        if(sqrt(σ2_max) > 0.099999 && sqrt(σ2_max) < 0.100001)
            println("Converge! A = ", A)
            break
        else
            if(global_ite ≤ 3)
                A = (0.1 / sqrt(σ2_max)) * A
            else
                A = (A + (0.1 / sqrt(σ2_max)) * A)*0.5
            end
            println("New A = ", A)

            # set force
            A, fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc = (
                Force_density.set_force_density(type_force, ρ_r, dρ_dr_r, r_c, r_o, R, A, j_0, j_1)
                )
            m_rr_co, m_rr_cr, m_rr_oc, m_rth_co, m_rth_cr, m_rth_oc = (
                Force_density.set_magnetic_stress(A, r_c, r_o, R, j_0, j_1)
                )
        end

    end

    σ2_max = Perturb_star.cal_σ!(σ, σr1, σr2, σr3, r_c, r_o, r_g, T1, T2, μ, ξr, ξt, nr)

    for i=1:nr
        if r_g[i] < r_c  # core
            fr_g[i] = fr_co(r_g[i])
            ft_g[i] = ft_co(r_g[i])
            M1[i] = m_rr_co(r_g[i])
            M2[i] = m_rth_co(r_g[i])
        elseif r_g[i] <= r_o # crust
            fr_g[i] = fr_cr(r_g[i])
            ft_g[i] = ft_cr(r_g[i])
            M1[i] = m_rr_cr(r_g[i])
            M2[i] = m_rth_cr(r_g[i])
        else # ocean 
            fr_g[i] = fr_oc(r_g[i])
            ft_g[i] = ft_oc(r_g[i])
            M1[i] = m_rr_oc(r_g[i])
            M2[i] = m_rth_oc(r_g[i])
        end
    end

    abs_eps = abs(ε_s - ε_f)
    σ_max = sqrt(σ2_max)

    out = open("result2.dat", "w")
    println(out, "#nr = $nr, M = $(M/M0) Mo,  R = $(R/1.0e5)[km], r_c = $r_c,  r_o = $r_o")
    println(out, "#ρc = $ρc, ρo = $ρo, ρ_min = $ρ_min") 
    println(out, "#Forcetype = $(type_force), A = $(A), j_0 = $j_0, j_1 = $j_1") 
    println(out, "#ε_f = $(ε_f), ε_s = $(ε_s), |ε_s - ε_f| = $(abs_eps), σ_max = $(σ2_max)") 



    @printf(out,   "%21s %21s %21s %21s %21s %21s %21s", "1(r_g[i])", "2(ρ[i])", "3(dρ_dr[i])", "4(cs2[i])", "5(μ[i])", "6(fr_g[i])", "7(ft_g[i])")
    @printf(out,   " %21s %21s %21s %21s %21s %21s %21s", "8(δρ_f[i])", "9(δϕ_f[i])", "10(δρ[i])", "11(δϕ[i])","12(M2[i])","13(ξr[i])","14(ξt[i])")
    @printf(out,   " %21s %21s %21s %21s %21s %21s \n", "15(T1[i])", "16(T2[i])", "17(σ[i])", "18(σr1[i])", "19(σr2[i])", "20(σr3[i])")


    for i=1:nr
        @printf(out, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e", r_g[i], ρ[i], dρ_dr[i], cs2[i], μ[i], fr_g[i], ft_g[i])
        @printf(out, " %.15e %.15e %.15e %.15e %.15e %.15e %.15e", δρ_f[i], δϕ_f[i], δρ[i], δϕ[i], M2[i], ξr[i], ξt[i])
        @printf(out, " %.15e %.15e %.15e %.15e %.15e %.15e \n", T1[i], T2[i], σ[i],  σr1[i], σr2[i], σr3[i])

        #        print(out, "$(r_g[i]) $(ρ[i]) $(fr_g[i]) $(ft_g[i])")
#        print(out, " $(δρ[i]) $(δp[i]) $(δϕ[i]) $(M1[i]) $(M2[i]) $(ξr[i]) $(ξt[i])")
#        print(out, " $(T1[i]) $(T2[i]) $(σ[i]) $(σr1[i]) $(σr2[i]) $(σr3[i]) \n")
    end
    close(out)
#    out = open("T1T2.dat","w")
#    for i=1:nr
#        print(out, "$(r_g[i])  $(δp[i] - T1[i]) $(δp[i])  $(T1[i])  $(T2[i])  $(δp[i] - T1[i] - M1[i]) $(T2[i] + M2[i]) $(M1[i])  $(M2[i])\n")
#    end
#    close(out)

    #plot by using matplotlib
    fig, ax = plt.subplots()
    ax.plot(r_g, T1)
    ax.grid(true)
    ax.set_xlabel("r [cm]")
    ax.set_xlim(r_c, radius)
    ax.set_ylabel("T1")
    plt.savefig("T1.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, δp - T1)
    ax.grid(true)
    ax.set_ylabel("δp - T1")
    ax.set_xlim(0.0, radius)
    plt.savefig("deltap_T1.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, T2)
    ax.grid(true)
    ax.set_ylabel("T2")
    ax.set_xlim(r_c, radius)
    plt.savefig("T2.pdf")
    plt.close()


    fig, ax = plt.subplots()
    ax.plot(r_g, M1)
    ax.grid(true)
    ax.set_ylabel("M1")
    ax.set_xlim(0.0, radius)
    plt.savefig("M1.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, M2)
    ax.grid(true)
    ax.set_ylabel("M2")
    ax.set_xlim(0.0, radius)
    plt.savefig("M2.pdf")
    plt.close()


    fig, ax = plt.subplots()
    ax.plot(r_g, δp - T1 - M1)
    ax.grid(true)
    ax.set_ylabel("δp - T1 - M1")
    ax.set_xlim(0.0, radius)
    plt.savefig("deltap_T1_2.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, T2 + M2)
    ax.grid(true)
    ax.set_ylabel("T2 + M_rth")
    ax.set_xlim(0.0, radius)
    plt.savefig("T2_2.pdf")
    plt.close()

    # M1 % M2

    fig, ax = plt.subplots()
    ax.plot(r_g, M1)
    ax.grid(true)
    ax.set_ylabel("M1")
    ax.set_xlim(0.0, radius)
    plt.savefig("M1.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, M2)
    ax.grid(true)
    ax.set_ylabel("M2")
    ax.set_xlim(0.0, radius)
    plt.savefig("M2.pdf")
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
    ax.set_xlim(r_min, radius)
    plt.savefig("delta_rho.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, δϕ)
    ax.grid(true)
    ax.set_ylabel("δϕ")
    ax.set_xlim(r_min, radius)
    plt.savefig("delta_phi.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(r_g, σ)
    ax.grid(true)
    ax.set_ylabel("|σ|")
    ax.set_xlim(9.98e5, radius)
    plt.savefig("sigma.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.plot(r_g, σr1)
    ax.plot(r_g, σr2)
    ax.plot(r_g, σr3)
    ax.grid(true)
    ax.set_ylabel("|σ|")
    ax.set_xlim(0.99r_o, radius)
    ax.set_ylim(1.0e-6, 0.1)
    plt.savefig("sigma123.pdf")
    plt.close()

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.plot(r_g, ρ)
    ax.grid(true)
    ax.set_ylabel("ρ")
    ax.set_xlim(0.999r_o, radius)
    plt.savefig("density.pdf")
    plt.close()

end
main()