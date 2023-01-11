using Plots
using DifferentialEquations
using Dierckx

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

    ρ0 = 2.0e15
    p0 = K*ρ0^Γ
    m0 = 0.0
    r_min = 0.1
    f0 = [m0, p0]
    rspan = (r_min, 2.0e6)
    prob0 = ODEProblem(back_ground, f0, rspan, param)
    cb     = DiscreteCallback(condition, affect!)
    bg_star = solve(prob0, Tsit5(), callback=cb, reltol = 1.0e-10, abstol = 1.0e-10)
#    bg_star = solve(prob0, Tsit5(), callback=cb)

    nr = size(bg_star.t, 1) - 1
    #radius
    R = bg_star.t[nr]
    M = bg_star.u[nr][1]
    r = zeros(Float64, nr)
    ρ = zeros(Float64, nr)
    p = zeros(Float64, nr)
    m = zeros(Float64, nr)
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

    ρ_r = Spline1D(r, ρ)

    #grid point
    ng = 100
    r_g = range(1, R*0.99999, ng)
    ρ_g = zeros(Float64, ng)
    dρ_g_dr = zeros(Float64, ng)
    for i=1:ng
        ρ_g[i] = ρ_r(r_g[i])
        dρ_g_dr[i] = derivative(ρ_r, r_g[i]; nu=1)
    end

    plot(r, dρ_dr)
    plot!(r_g, dρ_g_dr)

end 

main()
