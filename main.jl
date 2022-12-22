###calculate mountain on a neutron star
#set ρ_c = 1, G = 1

using Plots
using LinearAlgebra

function integ_simpson(S, Nr, dr)
    ans = 0.0
    for i = 1:2:Nr-2
        ans += dr*(S[i] + 4*S[i+1] + S[i+2]) / 3.0
    end
    return ans
end

function diff_r(f, dr, Nr)
    df_dr = zeros(Float64, Nr)
#    df_dr[1] = (-f[3] + 4*f[2] - 3*f[1]) / (2*dr)
    df_dr[1] = (f[2] - f[1]) / (dr)
    for i=2:Nr-1
        df_dr[i] = (f[i+1] - f[i-1]) / (2*dr)
    end
    df_dr[Nr] = (3*f[Nr] - 4*f[Nr-1] + f[Nr-2]) / (2*dr)
    return df_dr
end

function diff_r2(f, dr, Nr)
    df_dr2 = zeros(Float64, Nr)
    df_dr2[1] = (f[3]- 2*f[2] + f[1]) / dr^2
    for i=2:Nr-1
        df_dr2[i] = (f[i+1] - 2*f[i] + f[i-1]) / (dr^2)
    end
    df_dr2[Nr] = (f[Nr] - 2*f[Nr-1] + f[Nr-2]) / dr^2
    return df_dr2
end


function cal_ϕ(r, dr, ρ, Nr, ℓ, m)

    S   = zeros(Float64, Nr)
    phi = zeros(Float64, Nr)
    #gravitational potential
    for i=1:Nr
        for j=1:Nr
            if r[i] ≤ r[j]
               S[j] = 4π*r[i]^(ℓ) * r[j]^(-ℓ+1) * ρ[j]
            else
               S[j] = 4π*r[j]^(ℓ+2) * r[i]^(-ℓ-1) * ρ[j]
            end
        end
        phi[i] =-integ_simpson(S, Nr, dr)
    end
    if 0 < ℓ
        phi[1] = 0.0 
    end
    return phi
end

function main()

    #number of radial grid
    Nr = 129
    #background matter p = Kρ^Γ 
    ρ  = zeros(Float64, Nr)
    p  = zeros(Float64, Nr)
    cs2= zeros(Float64, Nr) #c_s^2 = K Γ ρ^{Γ - 1}
    ϕ  = zeros(Float64, Nr)
    K  = 1.0  #initial guess
    Γ  = 2.0
    #radial coordinate r
    r  = zeros(Float64, Nr)

    r_min = 0.0
    r_max = 1.0
    dr    = (r_max - r_min) / (Nr-1)

    #perturbations
    δρ = zeros(Float64, Nr)
    δϕ = zeros(Float64, Nr)

    #set background star
    #center
    r[1]   = r_min
    ρ[1]   = 1.0
    for i=2:Nr-1
        r[i] = r[i-1] + dr
        ρ[i] = sin(π*r[i]) / (π*r[i])  #Γ = 1 
    end
    #surface
    r[Nr] = r_max
    ρ[Nr] = 0.0  

    ϕ = cal_ϕ(r, dr, ρ, Nr, 0, 0)

 #   display(ϕ)

#    plot(r, ϕ, xlims=(0.0,1.0))
    K = (-ϕ[1] + ϕ[Nr]) / (Γ / (Γ - 1))

    for i=1:Nr
        p[i]   = K*ρ[i]^Γ
        cs2[i] = K*Γ*ρ[i]^(Γ-1)
    end
###end make background


###perturbations

    ft = zeros(Float64, Nr)

    for i=1:Nr
        ft[i] = -2*(1 - r[i])^2
    end
    Nstep = 1

    δ = zeros(Float64, Nr)
    for step = 1:Nstep

        for i=1:Nr
            δρ[i] = -(ρ[i]*δϕ[i] - ft[i] * r[i]) / cs2[i]
        end
        δρ[Nr] = 0.0
        δϕ = cal_ϕ(r, dr, δρ, Nr, 2, 2) / 2π
        for i=1:Nr
            δ[i] = δρ[i] * δϕ[i]
        end
    end

    dϕ_dr   = zeros(Float64, Nr)
    d2ϕ_dr2 = zeros(Float64, Nr)

    dϕ_dr   = diff_r(ϕ, dr, Nr)
    d2ϕ_dr2 = diff_r(dϕ_dr, dr, Nr)

    

    plot(r, d2ϕ_dr2)
end

main()

#A = [1.0 -2.0; -1.0 1.0]
#b = [1.0, -2.0]

#    display(A)
#    display(b)
#    c = A \ b
#    display(A*c)
#    display(inv(A)*A)
