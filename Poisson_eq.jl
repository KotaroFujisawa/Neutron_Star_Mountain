module Poisson_eq
using DifferentialEquations
using Dierckx
using QuadGK
    #solve perturb Poisson equation by shooting method
    function Poisson(δρ_r, ℓ, r_min, R, r, nr, max_ite=100000)
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
            probP = ODEProblem(Poisson, f0, rspan, param)
            ppo = solve(probP, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
        
            np = size(ppo.t, 1)
            #dδϕ/dr + (ℓ+1)/R δϕ = 0 at r = R 
            #-> f = -dδϕ/dr  /  (ℓ+1)/R * δϕ
            f2 = (-ppo.u[np][2] / ((ℓ+1) / (R) * ppo.u[np][1]) - 1.0)

            #Newton method
            fp = (f2 - f) / ϵ
            dδϕ_dr_c = dδϕ_dr_c - 1.0e3 * f / fp

            if(abs(f) < 1.0e-6) 
                print("Converge j = $j dδϕ_dr_c = $dδϕ_dr_c f = $f \n")
                for i = 1:nr
                    δϕ[i], dδϕ_dr[i] = ppo(r[i])
                end
                break
            end
        end
        return δϕ, dδϕ_dr
    end
end