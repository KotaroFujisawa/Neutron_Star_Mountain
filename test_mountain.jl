using Plots
using DifferentialEquations
using Dierckx
using QuadGK

include("background_star.jl")
using .background_star



function main()
    #simple NS model
    Γ = 2
    K = 46000     #cgs
    M0 = 1.989e33 #cgs
    G = 6.6743e-8 #cgs

    ℓ  = 2        #ℓ = m = 2 deformation
    β2 = ℓ*(ℓ+1)  #β^2 = ℓ(ℓ+1)

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
    nr, R, M, r, ρ, p, cs2, m, dρ_dr, dp_dr, d2ρ_dr2 = (
        background_star.make_bg_star_p(
            K, Γ, ρ0, ρ_min, r_min, r_max
        )
    )

    #background star
    print("grid number = $nr, M = $(M/M0)  R = $(R/1.0e5) [km] \n")

    #interpolation by spline curve
    ρ_r = Spline1D(r, ρ)
 
    r_c = background_star.cal_r_from_ρ(0.8*R,  ρc, ρ_r)
    r_o = background_star.cal_r_from_ρ(0.95*R, ρo, ρ_r)

    print("r_c = $r_c  \n")
    print("r_o = $r_o  \n")

    #perturbation
    δρ = zeros(Float64, nr)
    δϕ = zeros(Float64, nr)
    dδϕ_dr = zeros(Float64, nr)

### fluid star ###

    #force A f_i = -Aρ ∇_i (r^2 Y_lm)
    A = 1e9
    fr(r) = -2A*r*ρ_r(r)
    ft(r) =  -A*r*ρ_r(r)

    ###perturb Poisson eq
    rspan   = (r_min, R)
    param = ()
    δρ_r = Spline1D(r, δρ)

    function Poisson(du, u, param, r)
#        δϕ   = u[1]
#        δϕ_dr= u[2]
    
        β2 = 2*(2+1)
        du[1] = u[2]
        du[2] = -2/r*u[2] + β2/r^2*u[1] + 4π*G*δρ_r(r)
#        tmp = δρ_r(r)
#        print("r = $r, δρ = $tmp \n")
    end

    max_ite = 10
    #initial guess
    dδϕ_dr_c = -1.0e6
    for ite=1:max_ite
        
        for i=1:nr
            δρ[i] = -(ρ[i]*δϕ[i] - ft(r[i])) / cs2[i]
        end
        δρ_r = Spline1D(r, δρ)
        ϵ    = 1.0e-6
        S_int(r) = 4π*G*δρ_r(r) * 10r_min^(ℓ) * r^(-ℓ+1)
        S = quadgk(S_int, 10r_min, R)
        dδϕ_dr_c = -(S[1] ) / (10r_min)
#        print(dδϕ_dr_c, " \n")
        for j=1:10000
        #shooting
            f0 = [0.0, dδϕ_dr_c]

            probP = ODEProblem(Poisson, f0, rspan, param)
            ppo = solve(probP, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
        
            np = size(ppo.t, 1)
            #dδϕ/dr + (ℓ+1)/R δϕ = 0 at r = R 
            #-> f = -dδϕ/dr  /  (ℓ+1)/R * δϕ
            f = -ppo.u[np][2] / ((ℓ+1) / (R) * ppo.u[np][1])
            f -= 1.0
            f0 = [0.0, dδϕ_dr_c*(1.0 + ϵ)]
            #shooting
            probP = ODEProblem(Poisson, f0, rspan, param)
            ppo = solve(probP, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
        
            np = size(ppo.t, 1)
            #dδϕ/dr + (ℓ+1)/R δϕ = 0 at r = R 
            #-> f = -dδϕ/dr  /  (ℓ+1)/R * δϕ
            f2 = -ppo.u[np][2] / ((ℓ+1) / (R) * ppo.u[np][1])
            f2 -= 1.0
            #Newton method
            fp = (f2 - f) / ϵ
            dδϕ_dr_c = dδϕ_dr_c - 1.0e3 * f / fp

            if ite == 10000
#                print("j = $j dδϕ_dr_c = $dδϕ_dr_c f = $f \n")
            end
            if(abs(f) < 1.0e-6) 
#                print("j = $j dδϕ_dr_c = $dδϕ_dr_c f = $f \n")
                break
            end
        end
        f0 = [0.0, dδϕ_dr_c]
        probP = ODEProblem(Poisson, f0, rspan, param)
        ppo = solve(probP, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
        for i = 1:nr
            δϕ[i], dδϕ_dr[i] = ppo(r[i])
        end

    end

    #quadrupole 
    Q22(r) = δρ_r(r)*r^4

    I0 = 1.0e45

    S = quadgk(Q22, r_min, R)
    ε_f = sqrt(8π/15) * S[1] / I0
    print("ε_f = $ε_f \n")

    plot(r, ρ)

#solid crust

end

main()
