using Plots
using DifferentialEquations
using Dierckx
using QuadGK

function back_ground(du, u, param, r)
    G = 6.6743e-8
    Γ = param[1]
    K = param[2]
    #u[1] -> m(r), u[2] -> p(r)
    m = u[1]
    p = u[2]
    #p = Kρ^Γ
    ρ = (abs(p)/K)^(1/Γ) 
    du[1] = 4π * r^2 * ρ
    du[2] = -ρ * G * m / r^2
#    print("r = ", r, ", m =", m, ", ρ = ", p, "\n")
end

function condition(u, r, integrator)
    p_min = 1.0e9
    u[2] < p_min #integrate until u[2] = p = 0
end

function affect!(integrator)
    terminate!(integrator)
end

function main()
    #simple NS model
    Γ = 2
    K = 46000
    M0 = 1.989e33 #g
    G = 6.6743e-8
    param = (Γ, K)
    ℓ  = 2
    β2 = ℓ*(ℓ+1)

    ρ0 = 2.0e15   #central density     r = r_min
    ρc = 2.0e14   #core-crust boundary r = r_c
    ρo = 1.0e11   #crust-ocian boundary r = r_o 
    p0 = K*ρ0^Γ
    m0 = 0.0

    r_min = 0.1   #[cm]
    r_c = 0.0
    r_o = 0.0
    f0 = [m0, p0]
    rspan = (r_min, 2.0e6)
    prob0 = ODEProblem(back_ground, f0, rspan, param)
    cb     = DiscreteCallback(condition, affect!)
    bg_star = solve(prob0, Tsit5(), callback=cb, reltol = 1.0e-10, abstol = 1.0e-10)

    nr = size(bg_star.t, 1) - 1
    #radius
    R = bg_star.t[nr]
    M = bg_star.u[nr][1]
    r = zeros(Float64, nr)
    ρ = zeros(Float64, nr)
    p = zeros(Float64, nr)
    m = zeros(Float64, nr)

    δρ = zeros(Float64, nr)
    δϕ = zeros(Float64, nr)
    dδϕ_dr = zeros(Float64, nr)

    for i=1:nr
        r[i] = bg_star.t[i]
        p[i] = bg_star.u[i][2]
        ρ[i] = (p[i]/K)^(1/Γ)
        m[i] = bg_star.u[i][1]
    end 
    #background star
    print("n r = ", nr, " M = ", M/M0, " R = ", R/1.0e5, "\n")

    #dp_dr, dρ_dr
    dp_dr = zeros(Float64, nr)
    dρ_dr = zeros(Float64, nr)
    dρ_dr2= zeros(Float64, nr)
    cs2   = zeros(Float64, nr)
    r2    = zeros(Float64, nr)
    for i=1:nr
        dp_dr[i] = -ρ[i]*G*m[i]/r[i]^2
        dρ_dr[i] = -ρ[i]^(2-Γ) / (Γ*K) * G*m[i] /r[i]^2
        cs2[i]   = K*Γ*ρ[i]^(Γ-1)
        dρ_dr2[i]= (2ρ[i]^(2-Γ)/(Γ*K) * G*m[i]/r[i]^3
                    - (2-Γ)/(Γ*K) * ρ[i]^(1-Γ)*dρ_dr[i] * G*m[i]/r[i]^2
                    - ρ[i]^(3-Γ) / (Γ*K) * (4π*G)
                    )
    end

    #interpolation by spline curve
    ρ_r = Spline1D(r, ρ)
 
    r_c = 0.8e6
    for i=1:10
        r_c -= (ρ_r(r_c) - ρc) / derivative(ρ_r, r_c;nu=1)
    end
    r_o = 0.8e6
    for i=1:10
        r_o -= (ρ_r(r_o) - ρo) / derivative(ρ_r, r_o;nu=1)
    end

    print("r_c = $r_c  \n")
    print("r_o = $r_o  \n")

    #grid point
    ng = 1000
    r_g = range(10, 0.99999*R, ng)
    ρ_g = zeros(Float64, ng)
    dρ_g_dr = zeros(Float64, ng)
    dρ_g_dr2= zeros(Float64, ng)
    for i=1:ng
        ρ_g[i] = ρ_r(r_g[i])
        dρ_g_dr[i] = derivative(ρ_r, r_g[i]; nu=1)
        dρ_g_dr2[i]= derivative(ρ_r, r_g[i]; nu=2)
    end
###

### fluid star

    #force f_i = -Aρ ∇_i (r^2 Y_lm)
    A = 1e9
    fr(r) = -2A*r*ρ_r(r)
    ft(r) =  -A*r*ρ_r(r)

    ###perturbe Poisson eq
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
