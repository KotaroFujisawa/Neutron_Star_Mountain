module Poisson_eq
using DifferentialEquations
using Dierckx
using QuadGK
using BoundaryValueDiffEq

    #solve perturb Poisson equation by shooting method
    function Poisson(δρ_r, ℓ, r_min, R, r_g, nr, ϵ = 1.0e-6, max_ite=100000)
        G = 6.6743e-8 #cgs
        ϵ = 1.0e-6

        #ℓ = m = 2 Poisson
        function source(du, u, param, r)
#        δϕ   = u[1]
#        δϕ_dr= u[2]    
            β2 = ℓ*(ℓ+1)
            du[1] = u[2]
            du[2] = -2/r*u[2] + β2/r^2*u[1] + 4π*G*δρ_r(r)
        end

        #initial guess
        S_int(r) = 4π*G*δρ_r(r) * 10r_min^(ℓ) * r^(-ℓ+1)
        S = quadgk(S_int, 10r_min, R)

        dδϕ_dr_c = - (S[1]) / (10r_min)
        param = ()
        rspan = (r_min, R)

 #       println("dδϕ_dr_c = ", dδϕ_dr_c)
        δϕ = zeros(Float64, nr)
        dδϕ_dr = zeros(Float64, nr)

        for i=1:max_ite
            f0 = [0.0, dδϕ_dr_c]

            prob = ODEProblem(source, f0, rspan, param)
            ppo = solve(prob, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
        
            np = size(ppo.t, 1)
            #dδϕ/dr + (ℓ+1)/R δϕ = 0 at r = R 
            #-> f = -dδϕ/dr  /  (ℓ+1)/R * δϕ
            f = (-ppo.u[np][2] / ((ℓ+1) / (R) * ppo.u[np][1]) - 1.0)

            f0 = [0.0, dδϕ_dr_c*(1.0 + ϵ)]
            #shooting
            prob2 = ODEProblem(source, f0, rspan, param)
            ppo2 = solve(prob2, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
        
            np = size(ppo2.t, 1)
            #dδϕ/dr + (ℓ+1)/R δϕ = 0 at r = R 
            #-> f = -dδϕ/dr  /  (ℓ+1)/R * δϕ
            f2 = (-ppo2.u[np][2] / ((ℓ+1) / (R) * ppo2.u[np][1]) - 1.0)

            #Newton method
            fp = (f2 - f) / ϵ
            dδϕ_dr_c = dδϕ_dr_c - 5.0e1 * f / fp

            if(abs(f) < 1.0e-6) 
 #               print("Converge i = $i dδϕ_dr_c = $dδϕ_dr_c f = $f \n")
                for i = 1:nr
                    δϕ[i], dδϕ_dr[i] = ppo(r_g[i])
                end
                break
            end
        end
#        println("f:dδϕ_dr_c = ", dδϕ_dr_c)
        return δϕ, dδϕ_dr
    end

    function Poisson_BVP(δρ_r, ℓ, r_min, R, r_g, nr)
        G = 6.6743e-8 #cgs
        ϵ = 1.0e-6
    
        #ℓ = m = 2 Poisson
        function source(du, u, param, r)
    #        δϕ   = u[1]
    #        δϕ_dr= u[2]    
            β2 = ℓ*(ℓ+1)
            du[1] = u[2]
            du[2] = -2/r*u[2] + β2/r^2*u[1] + 4π*G*δρ_r(r)
        end
    
        function bc(res, sol, param, r)
            res[1] = sol(r_min)[1]
            res[2] = (-sol(R)[2] / ((ℓ+1) / (R) * sol(R)[1]) - 1.0)
        end

        #initial guess
        S_int(r) = 4π*G*δρ_r(r) * 10r_min^(ℓ) * r^(-ℓ+1)
        S = quadgk(S_int, 10r_min, R)
        
        dδϕ_dr_c = - (S[1]) / (10r_min)
        init = [0.0, dδϕ_dr_c]
        param = ()
        rspan = (r_min, R)
        bvp = BVProblem(source, bc, init, rspan)
        ppo = solve(bvp, Shooting(Vern6()))
        δϕ = zeros(Float64, nr)
        dδϕ_dr = zeros(Float64, nr)

        for i = 1:nr
            δϕ[i], dδϕ_dr[i] = ppo(r_g[i])
        end
        return δϕ, dδϕ_dr
    end
end