module Perturb_star
using Dierckx
using QuadGK

    function δρ_fluid_star(r_g, ρ, δϕ, cs2, fr, ft, nr)
        δρ = zeros(Float64, nr)
        for i=1:nr
            δρ[i] = -(ρ[i]*δϕ[i] - ft(r_g[i])) / cs2[i]
        end
        return δρ
    end

    function ellipticity(δρ_r, r_min, R, I0 = 1.0e45)
        Q22(r) = δρ_r(r) * r^4
    
        S = quadgk(Q22, r_min, R)
        ε = sqrt(8π/15) * S[1] / I0
        return ε
    end

    function cal_ξ_T(ρ, dρ_dr, d2ρ_dr2, cs2, dcs2_dr, μ, dμ_dr, δϕ, dδϕ_dr, β2)

        β = sqrt(β2)

        function source(du, u, param, r)
            #u[1] = ξr, u[2] = ξt, u[3] = T1, u[4] = T2
            ξr = u[1]
            ξt = u[2]
            U1 = u[3]
            U2 = u[4]

            du[1] = ξr/r - β/(2r)*ξt + 3/(4μ(r))*T1
            du[2] =-β/r*ξr + ξt/r + β/μ(r)*T2
            #denominator
            denom = (1 + 3cs2(r)*ρ(r) / (4μ(r)))
            du[3] = (
                ρ(r)*δϕ(r) - fr(r) 
                -(dcs2_dr(r)*(3ρ(r)+r*dρ_dr(r)) + cs2(r)*(3β2*ρ(r)/(2r) + dρ_dr(r) - r*dρ_dr(r)^2/ρ + r*d2ρ_dr2(r)))*ξr/r
                +(dcs2_dr(r)*3ρ(r) + cs2(r)*(3ρ(r)/r + dρ_dr(r)))*β/(2r)*ξt
                -(3/r + dcs2_dr(r)*3ρ(r)/(4μ(r))+cs2(r)*(3ρ(r)/r - ρ(r)*dμ_dr(r)/μ(r) + dρ_dr(r))*3/(4μ(r)))*T1
                +(1 +3cs2(r)*ρ(r) / (2μ(r)))*beta2/r * T2
             ) / denom
            du[4] = (
                ρ(r) / r * δϕ(r) - ft(r) 
                - cs2(r)*(3ρ(r) + r*dρ_dr(r))*ξr/r^2
                + (3cs2(r)*ρ(r) / 2 + (1-2/β2)*μ(r))*β/r^2*ξt
                + (1/2 - 3cs2(r)*ρ/(4μ(r)))*T1/r - 3T2/r
            )
        end
        

    end

end
