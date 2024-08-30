module Poisson_eq
using DifferentialEquations
using Dierckx
using QuadGK
using BoundaryValueDiffEq
using ForwardDiff

    function Poisson_BVP(δρ_r, ℓ, r_min, R, r_g, nr, guess)
        G = 6.6743e-8 #cgs
        ϵ = 1.0e-6
    
        
        #ℓ = m = 2 Poisson
        function source(du, u, param, r)
    #        δϕ   = u[1]
    #        δϕ_dr= u[2]    
            β2 = ℓ*(ℓ+1)
            du[1] = u[2]
            du[2] = -2/r*u[2] + β2/r^2*u[1] + 4π*G*δρ_r(r)
            nothing
        end
    
        function bc(res, sol, param, r)
            res[1] = sol(r_min)[1]
            res[2] = (-sol(R)[2] / ((ℓ+1) / (R) * sol(R)[1]) - 1.0)
            nothing
        end

        #initial guess
        S_int(r) = 4π*G*δρ_r(r) * 10r_min^(ℓ) * r^(-ℓ+1)
        S = quadgk(S_int, 10r_min, R)
        dδϕ_dr_c = - (S[1]) / (10r_min)
#        print(dδϕ_dr_c, " ", guess, "\n")
        if(abs(guess) < 1.0)
            init = [0.0, dδϕ_dr_c]
        else
            init = [0.0, guess]
        end
        param = ()
        rspan = (r_min, R)
        bvp = BVProblem(source, bc, init, rspan, reltol = 1.0e-10, abstol = 1.0e-10)

#        ppo = solve(bvp, Shooting(Vern6()), save_everystep = false)
        ppo = solve(bvp, Shooting(Vern6()))
        δϕ = zeros(Float64, nr)
        dδϕ_dr = zeros(Float64, nr)

        for i = 1:nr
            δϕ[i], dδϕ_dr[i] = ppo(r_g[i])
        end
        return δϕ, dδϕ_dr
    end
end