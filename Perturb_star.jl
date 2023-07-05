module Perturb_star
using Dierckx
using QuadGK
using DifferentialEquations
using BoundaryValueDiffEq

    function δρ_fluid_star(r_g, ρ, δϕ, cs2, fr, ft, nr)
        δρ = zeros(Float64, nr)
        for i=1:nr
            δρ[i] = -(ρ[i]*δϕ[i] - ft(r_g[i])*r_g[i]) / cs2[i]
        end
        return δρ
    end

    function ellipticity(δρ_r, r_min, R, I0 = 1.0e45)
        Q22(r) = δρ_r(r) * r^4
    
        S = quadgk(Q22, r_min, R)
        ε = sqrt(8π/15) * S[1] / I0
        return ε
    end

    function cal_ξ_T_BVP(ρ, dρ_dr, d2ρ_dr2, cs2, dcs2_dr, 
        μ, dμ_dr, δϕ, dδϕ_dr, fr, ft, β2, r_c, r_o, r_g, nr)

        β = sqrt(β2)
        println("BVP")
        #four ODE for T1, T2, ξr, ξt
        #BC 1&2, T2 = 0 at r_c & r_o
        #BC 3&4 The relation T1, ξr, ξt at r_c & r_o
        function source(du, u, param, r)
            #u[1] = ξr, u[2] = ξt, u[3] = T1, u[4] = T2
            ξr = u[1]
            ξt = u[2]
            T1 = u[3]
            T2 = u[4]
            du[1] = ξr/r - β/(2r)*ξt + 3/(4μ(r))*T1
            du[2] =-β/r*ξr + ξt/r + β/μ(r)*T2
            #denominator
            denom = (1 + 3cs2(r)*ρ(r) / (4μ(r)))
            du[3] = (
                ρ(r)*dδϕ_dr(r) - fr(r) 
                -(dcs2_dr(r)*(3ρ(r)+r*dρ_dr(r)) + cs2(r)*(3β2*ρ(r)/(2r) + dρ_dr(r) - r*dρ_dr(r)^2/ρ(r) + r*d2ρ_dr2(r)))*ξr/r
                +(dcs2_dr(r)*3ρ(r) + cs2(r)*(3ρ(r)/r + dρ_dr(r)))*β/(2r)*ξt
                -(3/r + dcs2_dr(r)*3ρ(r)/(4μ(r))+cs2(r)*(3ρ(r)/r - ρ(r)*dμ_dr(r)/μ(r) + dρ_dr(r))*3/(4μ(r)))*T1
                +(1 +3cs2(r)*ρ(r) / (2μ(r)))*β2/r * T2
             ) / denom
            du[4] = (
                ρ(r) / r * δϕ(r) - ft(r) 
                - cs2(r)*(3ρ(r) + r*dρ_dr(r))*ξr/r^2
                + (3cs2(r)*ρ(r) / 2 + (1-2/β2)*μ(r))*β/r^2*ξt
                + (1/2 - 3cs2(r)*ρ(r)/(4μ(r)))*T1/r - 3T2/r
            )
        end
            
        function bc(res, sol, param, r)
            ξra = sol(r_c)[1]
            ξta = sol(r_c)[2]
            T1a = ((ρ(r_c)*δϕ(r_c) - ft(r_c)*r_c - 
                cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξra - 3β*ρ(r_c)/2r_c*ξta))
                / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
            T2a = sol(r_c)[4]

            ξrb = sol(r_o)[1]
            ξtb = sol(r_o)[2]
            T1b = ((ρ(r_o)*δϕ(r_o) - ft(r_o)*r_o - 
                cs2(r_o)*( (3ρ(r_o)/r_o + dρ_dr(r_o))*ξrb - 3β*ρ(r_o)/2r_o*ξtb))
                / (1 + 3cs2(r_o)*ρ(r_o)/(4μ(r_o))))
            T2b = sol(r_o)[4]
            
            res[1] = T2a 
            res[2] = (T1a - sol(r_c)[3]) 
            res[3] = T2b 
            res[4] = (T1b - sol(r_o)[3])
        end

        rspan = (r_c, r_o)
        #initial guess at r_c
        ξr_0 = -400.0
        ξt_0 =  -10.0
        T2_0 = 0.0
        T1_0 = ((ρ(r_c)*δϕ(r_c) - ft(r_c)*r_c - 
            cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξr_0 - 3β*ρ(r_c)/2r_c*ξt_0))
            / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
 
        init =[ξr_0, ξt_0, T1_0, T2_0]
        println(init)
        bvp  = BVProblem(source, bc, init, rspan)
        ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-8, abstol = 1.0e-8)

        ξr_1 = ans(r_c)[1]
        ξt_1 = ans(r_c)[2]
        init_ans = [ξr_1, ξt_1]    

        println(init_ans)

        T1_g = zeros(Float64, nr)
        T2_g = zeros(Float64, nr)
        ξr_g = zeros(Float64, nr)
        ξt_g = zeros(Float64, nr)

   #     display(ans)
        for i = 1:nr
            if(r_c ≤ r_g[i] && r_g[i] ≤ r_o)
                ξr_g[i], ξt_g[i], T1_g[i], T2_g[i] = ans(r_g[i])
            else
                ξr_g[i], ξt_g[i], T1_g[i], T2_g[i] = 0.0, 0.0, 0.0, 0.0
            end
        end
        return T1_g, T2_g, ξr_g, ξt_g
    end
    
end