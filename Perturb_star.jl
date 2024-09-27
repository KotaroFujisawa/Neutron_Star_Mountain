module Perturb_star
using NLsolve
using Dierckx
using QuadGK
using DifferentialEquations
using BoundaryValueDiffEq
using ForwardDiff

    function δρ_fluid_star(r_c, r_o, r_g, ρ, δϕ, cs2, 
        fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc, nr)
        δρ = zeros(Float64, nr)
        for i=1:nr
            if r_g[i] < r_c  # core
                δρ[i] = -(ρ[i]*δϕ[i] - ft_co(r_g[i])*r_g[i]) / cs2[i]
            elseif r_g[i] <= r_o # crust
                δρ[i] = -(ρ[i]*δϕ[i] - ft_cr(r_g[i])*r_g[i]) / cs2[i]
            else # ocean
                δρ[i] = -(ρ[i]*δϕ[i] - ft_oc(r_g[i])*r_g[i]) / cs2[i]
            end
        end
        return δρ
    end

    function δρ_solid_star(r_c, r_o, r_g, ρ, dρ_dr, δϕ, cs2, 
        fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc, ξr, ξt, T1, μ, β, nr)
        δρ = zeros(Float64, nr)
        for i=1:nr
            if r_g[i] < r_c # core
                δρ[i] =  -(ρ[i]*δϕ[i] - ft_co(r_g[i])*r_g[i]) / cs2[i]
            elseif r_g[i] <= r_o # crust
                δρ[i] = -(3ρ[i]/r_g[i] + dρ_dr[i])*ξr[i] + 3β*ρ[i]/(2r_g[i])*ξt[i] - 3ρ[i]/(4μ[i])*T1[i]
            else # ocean
                δρ[i] =  -(ρ[i]*δϕ[i] - ft_oc(r_g[i])*r_g[i]) / cs2[i]
            end
        end
        return δρ
    end

    function ellipticity(δρ_r, r_min, R, I0 = 1.0e45)
        Q22(r) = δρ_r(r) * r^4
    
        S = quadgk(Q22, r_min, R)
        ε = sqrt(8π/15) * S[1] / I0
        return ε
    end

    function cal_ξ_T_BVP0!(ρ, dρ_dr, d2ρ_dr2, cs2, dcs2_dr, 
        μ, dμ_dr, δϕ, dδϕ_dr, fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc,
        m_rr_co, m_rr_cr, m_rr_oc, m_rth_co, m_rth_cr, m_rth_oc,
        β2, r_c, r_o, r_g, nr, guess2, type_BC)

        β = sqrt(β2)
        println("BVP0")

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
                ρ(r)*dδϕ_dr(r) - fr_cr(r) 
                -(dcs2_dr(r)*(3ρ(r)+r*dρ_dr(r)) + cs2(r)*(3β2*ρ(r)/(2r) + dρ_dr(r) - r*dρ_dr(r)^2/ρ(r) + r*d2ρ_dr2(r)))*ξr/r
                +(dcs2_dr(r)*3ρ(r) + cs2(r)*(3ρ(r)/r + dρ_dr(r)))*β/(2r)*ξt
                -(3/r + dcs2_dr(r)*3ρ(r)/(4μ(r))+cs2(r)*(3ρ(r)/r - ρ(r)*dμ_dr(r)/μ(r) + dρ_dr(r))*3/(4μ(r)))*T1
                +(1 +3cs2(r)*ρ(r) / (2μ(r)))*β2/r * T2
             ) / denom
            du[4] = (
                ρ(r) / r * δϕ(r) - ft_cr(r) 
                - cs2(r)*(3ρ(r) + r*dρ_dr(r))*ξr/r^2
                + (3cs2(r)*ρ(r) / 2 + (1-2/β2)*μ(r))*β/r^2*ξt
                + (1/2 - 3cs2(r)*ρ(r)/(4μ(r)))*T1/r - 3T2/r
            )
            nothing
        end
        
        # a: core-crsut boundary b: crust-ocean boundary

        # normal boundary condition
        function bc(res, sol, param, r)
            ξra = sol(r_c)[1]
            ξta = sol(r_c)[2]
            T1a = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
                cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξra - 3β*ρ(r_c)/2r_c*ξta))
                / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
            T2a = sol(r_c)[4]

            ξrb = sol(r_o)[1]
            ξtb = sol(r_o)[2]
            T1b = ((ρ(r_o)*δϕ(r_o) - ft_oc(r_o)*r_o - 
                cs2(r_o)*( (3ρ(r_o)/r_o + dρ_dr(r_o))*ξrb - 3β*ρ(r_o)/2r_o*ξtb))
                / (1 + 3cs2(r_o)*ρ(r_o)/(4μ(r_o))))
            T2b = sol(r_o)[4]
            
            res[1] = T2a 
            res[2] = (T1a - sol(r_c)[3]) 
            res[3] = T2b 
            res[4] = (T1b - sol(r_o)[3])
            nothing
#            println(res[3])
#            println("res = $res")
        end

        rspan = (r_c, r_o)
        function source2(x)
            ξr_0 = x[1]
            ξt_0 = x[2]

            T1_0 = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
            cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξr_0 - 3β*ρ(r_c)/2r_c*ξt_0))
            / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
            T2_0 = 0.0

            prob = ODEProblem(source, [ξr_0, ξt_0, T1_0, T2_0], rspan)
            sol = solve(prob, Vern7(), reltol = 1.0e-10, abstol = 1.0e-10)

            ξrb = sol(r_o)[1]
            ξtb = sol(r_o)[2]
            T1b = ((ρ(r_o)*δϕ(r_o) - ft_oc(r_o)*r_o - 
                cs2(r_o)*( (3ρ(r_o)/r_o + dρ_dr(r_o))*ξrb - 3β*ρ(r_o)/2r_o*ξtb))
                / (1 + 3cs2(r_o)*ρ(r_o)/(4μ(r_o))))
            T2b = sol(r_o)[4]
            
            T2c = abs(sol( (r_o + r_c)*0.5 )[4])
            T1c = abs(sol( (r_o + r_c)*0.5 )[3])

            res1 = T2b 
            res2 = (T1b - sol(r_o)[3]) 
            return [res1, res2]
        end

        ξr0 = -10000.0
        ξt0 = -300000.0

        #initial guess at r_c
        if(abs(guess2[1]) > 0.0)
            ξr0 = guess2[1]
            ξt0 = guess2[2]
        end
        T1_0 = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
        cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξr0 - 3β*ρ(r_c)/2r_c*ξt0))
        / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
        T2_0 = 0.0

        init = [ξr0, ξt0, T1_0, T2_0] 
        println(init)

        # shooting
        nsol = nlsolve(source2, [ξr0, ξt0], ftol=1.0e-6, autodiff = :forward, iterations = 1000)

        # get initial values at core-crust boundary
        ξr_0 = nsol.zero[1]
        ξt_0 = nsol.zero[2]

        T1_0 = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
        cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξr_0 - 3β*ρ(r_c)/2r_c*ξt_0))
        / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
        T2_0 = 0.0

        init = [ξr_0, ξt_0, T1_0, T2_0] 

        prob = ODEProblem(source, init, rspan)
        ans = solve(prob, Vern7(), reltol = 1.0e-10, abstol = 1.0e-10)

        ξr_1 = ans(r_c)[1]
        ξt_1 = ans(r_c)[2]
        T1_1 = ans(r_c)[3]
        T2_2 = ans(r_c)[4]
        init_ans = [ξr_1, ξt_1, T1_1, T2_2]    
        init2 = [ξr_1, ξt_1, T1_1]
        println(init_ans)

        T1_g = zeros(Float64, nr)
        T2_g = zeros(Float64, nr)
        ξr_g = zeros(Float64, nr)
        ξt_g = zeros(Float64, nr)

        for i = 1:nr
            if(r_c ≤ r_g[i] && r_g[i] ≤ r_o)
                ξr_g[i], ξt_g[i], T1_g[i], T2_g[i] = ans(r_g[i])
            else
                ξr_g[i], ξt_g[i], T1_g[i], T2_g[i] = 0.0, 0.0, 0.0, 0.0
            end
        end

        return T1_g, T2_g, ξr_g, ξt_g, init2
    end

    function cal_ξ_T_BVP!(ρ, dρ_dr, d2ρ_dr2, cs2, dcs2_dr, 
        μ, dμ_dr, δϕ, dδϕ_dr, fr_co, fr_cr, fr_oc, ft_co, ft_cr, ft_oc,
        m_rr_co, m_rr_cr, m_rr_oc, m_rth_co, m_rth_cr, m_rth_oc,
        β2, r_c, r_o, r_g, nr, guess2, type_BC)

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
                ρ(r)*dδϕ_dr(r) - fr_cr(r) 
                -(dcs2_dr(r)*(3ρ(r)+r*dρ_dr(r)) + cs2(r)*(3β2*ρ(r)/(2r) + dρ_dr(r) - r*dρ_dr(r)^2/ρ(r) + r*d2ρ_dr2(r)))*ξr/r
                +(dcs2_dr(r)*3ρ(r) + cs2(r)*(3ρ(r)/r + dρ_dr(r)))*β/(2r)*ξt
                -(3/r + dcs2_dr(r)*3ρ(r)/(4μ(r))+cs2(r)*(3ρ(r)/r - ρ(r)*dμ_dr(r)/μ(r) + dρ_dr(r))*3/(4μ(r)))*T1
                +(1 +3cs2(r)*ρ(r) / (2μ(r)))*β2/r * T2
             ) / denom
            du[4] = (
                ρ(r) / r * δϕ(r) - ft_cr(r) 
                - cs2(r)*(3ρ(r) + r*dρ_dr(r))*ξr/r^2
                + (3cs2(r)*ρ(r) / 2 + (1-2/β2)*μ(r))*β/r^2*ξt
                + (1/2 - 3cs2(r)*ρ(r)/(4μ(r)))*T1/r - 3T2/r
            )
            nothing
        end
            
        # a: core-crsut boundary b: crust-ocean boundary
        # normal boundary condition
        function bc(res, sol, param, r)
            ξra = sol(r_c)[1]
            ξta = sol(r_c)[2]
            T1a = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
                cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξra - 3β*ρ(r_c)/2r_c*ξta))
                / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
            T2a = sol(r_c)[4]

            ξrb = sol(r_o)[1]
            ξtb = sol(r_o)[2]
            T1b = ((ρ(r_o)*δϕ(r_o) - ft_oc(r_o)*r_o - 
                cs2(r_o)*( (3ρ(r_o)/r_o + dρ_dr(r_o))*ξrb - 3β*ρ(r_o)/2r_o*ξtb))
                / (1 + 3cs2(r_o)*ρ(r_o)/(4μ(r_o))))
            T2b = sol(r_o)[4]
            
            res[1] = T2a 
            res[2] = (T1a - sol(r_c)[3]) 
            res[3] = T2b 
            res[4] = (T1b - sol(r_o)[3])
            nothing
#            println(res[3])
#            println("res = $res")
        end

        # δp - T1 = 0 at core-crust boundary
        function bc2(res, sol, param, r)
            ξra = sol(r_c)[1]
            ξta = sol(r_c)[2]
            T1a = ((
                -cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξra - 3β*ρ(r_c)/2r_c*ξta))
                / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
            T2a = sol(r_c)[4]

            ξrb = sol(r_o)[1]
            ξtb = sol(r_o)[2]
            T1b = ((ρ(r_o)*δϕ(r_o) - ft_oc(r_o)*r_o - 
                cs2(r_o)*( (3ρ(r_o)/r_o + dρ_dr(r_o))*ξrb - 3β*ρ(r_o)/2r_o*ξtb))
                / (1 + 3cs2(r_o)*ρ(r_o)/(4μ(r_o))))
            T2b = sol(r_o)[4]
            
            res[1] = T2a 
            res[2] = (T1a - sol(r_c)[3])
            res[3] = T2b 
            res[4] = (T1b - sol(r_o)[3])
            nothing
        end

        function bc3(res, sol, param, r)
            ξra = sol(r_c)[1]
            ξta = sol(r_c)[2]
            T1a = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
                cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξra - 3β*ρ(r_c)/2r_c*ξta)
                + (m_rr_co(r_c) -  m_rr_cr(r_c)))
                / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
            T2a = m_rth_co(r_c) - m_rth_cr(r_c)

            ξrb = sol(r_o)[1]
            ξtb = sol(r_o)[2]
            T1b = ((ρ(r_o)*δϕ(r_o) - ft_oc(r_o)*r_o - 
                cs2(r_o)*( (3ρ(r_o)/r_o + dρ_dr(r_o))*ξrb - 3β*ρ(r_o)/2r_o*ξtb)
                + ( m_rr_oc(r_o) - m_rr_cr(r_o)))
                / (1 + 3cs2(r_o)*ρ(r_o)/(4μ(r_o))))
            T2b = m_rth_oc(r_o) - m_rth_cr(r_o)

            res[1] = T2a - sol(r_c)[4]
            res[2] = (T1a - sol(r_c)[3]) 
            res[3] = T2b - sol(r_o)[4]
            res[4] = (T1b - sol(r_o)[3])
            nothing
#            println("res = $res")
        end

        rspan = (r_c, r_o)

        #initial guess at r_c
        ξr_0 = -600.0
        ξt_0 =  50.0
        T2_0 = m_rth_co(r_c) - m_rth_cr(r_c)

        ξr_0 =  900.0
        ξt_0 = -6000.0

        ξr_0 = -700
        ξt_0 =  7000

#        ξr_0 = -25000
#        ξt_0 = 463664

        T1_0 = ((ρ(r_c)*δϕ(r_c) - ft_co(r_c)*r_c - 
            cs2(r_c)*( (3ρ(r_c)/r_c + dρ_dr(r_c))*ξr_0 - 3β*ρ(r_c)/2r_c*ξt_0))
            / (1 + 3cs2(r_c)*ρ(r_c)/(4μ(r_c))))
        if(abs(guess2[1]) > 0.0)
            ξr_0 = guess2[1]
            ξt_0 = guess2[2]
            T1_0 = guess2[3]
        end
        init =[ξr_0, ξt_0, T1_0, T2_0]
        println(init)
        if(type_BC==1) 
            bvp  = BVProblem(source, bc, init, rspan)
        end
        if(type_BC==2) 
            bvp  = BVProblem(source, bc2, init, rspan)
        end
        if(type_BC==3)
            bvp  = BVProblem(source, bc3, init, rspan)
        end
#        ans  = solve(bvp, Shooting(Vern7()))

    #    ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-5, abstol = 1.0e-5) 
        ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-6, abstol = 1.0e-6)
#        ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-7, abstol = 1.0e-7)
#        ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-7, abstol = 1.0e-7, save_everystep = false)
#        ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-8, abstol = 1.0e-8)
#        ans  = solve(bvp, Shooting(Vern6()), reltol = 1.0e-9, abstol = 1.0e-9)

        ξr_1 = ans(r_c)[1]
        ξt_1 = ans(r_c)[2]
        T1_1 = ans(r_c)[3]
        T2_2 = ans(r_c)[4]
        init_ans = [ξr_1, ξt_1, T1_1, T2_2]    
        init2 = [ξr_1, ξt_1, T1_1]
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

        return T1_g, T2_g, ξr_g, ξt_g, init2
    end
    
    
    function cal_σ!(σ, σr1, σr2, σr3, r_c, r_o, r_g, T1, T2, μ, ξr, ξt, nr)

        nθ = 100
        nφ = 100
    
        σ2  = zeros(Float64, nr, nθ, nφ)

        θ = range(0.0, stop=π,  length=nθ)
        φ = range(0.0, stop=2π, length=nφ)

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

        σ2_max = maximum(σ2)
        #maximum strain
        println("σ_max = ", sqrt(σ2_max))
    
        for i=1:nr
            σ[i] = sqrt(maximum(σ2[i,:,:]))
        end

        return σ2_max

    end
end