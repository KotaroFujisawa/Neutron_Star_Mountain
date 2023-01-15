module background_star
using DifferentialEquations
using Dierckx
#make polytropic star
    function make_bg_star_p(K, Γ, ρ0, ρ_min, r_min, r_max)
        G = 6.6743e-8
        p_min = K*ρ_min^Γ

        function condition(u, r, integrator)
            u[2] < p_min #integrate if u[2] = p > p_min
        end
        
        function affect!(integrator)
            terminate!(integrator)
        end
        
        function euler(du, u, param, r)
            #u[1] -> m(r), u[2] -> p(r)
            m = u[1]
            p = u[2]
            #p = Kρ^Γ
            ρ = (abs(p)/K)^(1/Γ) 
            du[1] = 4π * r^2 * ρ
            du[2] = -ρ * G * m / r^2
        end

        p0 = K*ρ0^Γ
        m0 = 0.0
        u0 = [m0, p0]
        param = [0.0, 0.0]
        rspan = (r_min, r_max)
        prob = DifferentialEquations.ODEProblem(euler, u0, rspan, param)
        cb     = DiscreteCallback(condition, affect!)
        bg_star = solve(prob, Tsit5(), callback=cb, reltol = 1.0e-10, abstol = 1.0e-10)
        #number of radial grids
        nr = size(bg_star.t, 1) - 1

        #make grid derivative
        R = bg_star.t[nr]
        M = bg_star.u[nr][1]

        r = zeros(Float64, nr)
        ρ = zeros(Float64, nr)
        p = zeros(Float64, nr)
        cs2= zeros(Float64,nr)
        m = zeros(Float64, nr)
        #dp_dr, dρ_dr
        dp_dr  = zeros(Float64, nr)
        dρ_dr  = zeros(Float64, nr)
        d2ρ_dr2= zeros(Float64, nr)


        for i=1:nr
            r[i] = bg_star.t[i]
            p[i] = bg_star.u[i][2]
            ρ[i] = (p[i]/K)^(1/Γ)
            cs2[i] = K*Γ*ρ[i]^(Γ-1)
            m[i] = bg_star.u[i][1]
            dp_dr[i] = -ρ[i]*G*m[i]/r[i]^2
            dρ_dr[i] = -ρ[i]^(2-Γ) / (Γ*K) * G*m[i] /r[i]^2
            d2ρ_dr2[i]= (2ρ[i]^(2-Γ)/(Γ*K) * G*m[i]/r[i]^3
                        - (2-Γ)/(Γ*K) * ρ[i]^(1-Γ)*dρ_dr[i] * G*m[i]/r[i]^2
                        - ρ[i]^(3-Γ) / (Γ*K) * (4π*G)
                        )
       
        end
        return nr, R, M, r, ρ, p, cs2, m, dρ_dr, dp_dr, d2ρ_dr2
    end
    
    function cal_r_from_ρ(r_i, ρc, ρ_r)

        max_ite = 50
        ϵ = 1.0e-8
        #initial guess
        r_c = r_i
        for i=1:max_ite
            r_c -= (ρ_r(r_c) - ρc) / derivative(ρ_r, r_c)
            if abs(ρ_r(r_c) - ρc) / ρc < ϵ
                break
            end
            if i==max_ite
                print("Not converge \n")
            end
        end
        return r_c
    end
end