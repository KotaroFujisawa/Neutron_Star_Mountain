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
    return phi / (2ℓ+1)
end

function main()

    #number of radial grid
    Nr = 129
    #number of rdial grid for coordinate
    Ncr = 12
    Nco = Nr - Ncr


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

    β = sqrt(2*(2+1))

    #perturbations
    δρ = zeros(Float64, Nr)
    δϕ = zeros(Float64, Nr)

    #shear modulas
    μ  = zeros(Float64, Nr)

    #stress
    ξr = zeros(Float64, Nr)
    ξt = zeros(Float64, Nr)
    T1 = zeros(Float64, Nr)
    T2 = zeros(Float64, Nr)

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

###crust###

κ = 1.0e-3
for i=Nco:Nr
    μ[i] =κ*ρ[i]
end
    

###perturbations

    ft = zeros(Float64, Nr)
    fr = zeros(Float64, Nr)
### external force ###
    for i=1:Nr
        ft[i] = -0.01*(1 - r[i])^2
        fr[i] =  0.01*ρ[i]*r[i]
    end
######################

###iteration fluid star### 
    Nstep = 30

    for step = 1:Nstep

        for i=1:Nr
            δρ[i] = -(ρ[i]*δϕ[i] - ft[i] * r[i]) / cs2[i]
        end
        δρ[Nr] = 0.0
        δϕ = cal_ϕ(r, dr, δρ, Nr, 2, 2)
    end

###iteration solid star

    dϕ_dr   = zeros(Float64, Nr)
    d2ϕ_dr2 = zeros(Float64, Nr)

    dδϕ_dr  = zeros(Float64, Nr)
    d2δϕ_dr2= zeros(Float64, Nr)

    dρ_dr   = zeros(Float64, Nr)
    d2ρ_dr2 = zeros(Float64, Nr)
    dcs2_dr = zeros(Float64, Nr)

    dμ_dr   = zeros(Float64, Nr)

    for step = 1:1

        ξr[Nco] = 0.0
        ξt[Nco] = 0.0
        T1[Nco] = 0.0
        T2[Nco] = 0.0

        dρ_dr = diff_r(ρ, dr, Nr)
        d2ρ_dr2 = diff_r(dρ_dr, dr, Nr)
        dμ_dr   = diff_r(μ, dr, Nr)
        for i=Nco:Nr-2
            ξr[i+1] = (
                ξr[i] 
                + dr * (ξr[i]/r[i] - β/(2r[i])*ξt[i] + 3/(4μ[i])*T1[i] )
            )
            ξt[i+1] = (
                ξt[i] 
                + dr * (-β/r[i]*ξr[i] + ξt[i]/r[i] + β/(μ[i])*T2[i] )
            )
            T1[i+1] = (
                T1[i] + dr / (1 + 3cs2[i]*ρ[i]/(4μ[i])) * 
                ( 
                    ρ[i] * dδϕ_dr[i] - fr[i] - 
                    ( 
                        + dcs2_dr[i] * (3ρ[i] + r[i]*dρ_dr[i])
                        + cs2[i] * 
                        (
                            3β^2*ρ[i]/(2*r[i]) + dρ_dr[i] 
                            -r[i]*dρ_dr[i]^2/ρ[i] + d2ρ_dr2[i] 
                        )
                    )*β/(2r[i])*ξt[i] +
                    (
                        dcs2_dr[i]^2*3ρ[i] + cs2[i]*
                        (
                            3ρ[i]/r[i] + dρ_dr[i]
                        )
                    )*β/(2r[i])*ξt[i] -
                    (
                        3/r[i] + dcs2_dr[i]*3ρ[i]/(4μ[i])
                        + cs2[i]*
                        (
                            3ρ[i]/r[i] - ρ[i]*dμ_dr[i]/μ[i] + dρ_dr[i]
                        ) * 3/(4μ[i])
                    ) * T1[i] +
                    (
                        1 + 3cs2[i]*ρ[i] / (2μ[i])
                    ) * β^2 / r[i] * T2[i]
                )
            )
            T2[i+1] = (
                T2[i] + dr * (
                    ρ[i]/r[i]*δϕ[i] - ft[i] 
                    - cs2[i]*(3ρ[i]+r[i]*dρ_dr[i])/r[i]^2*ξr[i]
                    + (3cs2[i]*ρ[i]/2 + (1 - 2/β^2)*μ[i])*β/r[i]^2*ξr[i]
                    + (1/2 - 3cs2[i]*ρ[i]/(4μ[i]))/r[i]*T1[i] 
                    - 3T2[i]/r[i]
                )
            )
        end
    end

    dϕ_dr   = diff_r(ϕ, dr, Nr)
    d2ϕ_dr2 = diff_r(dϕ_dr, dr, Nr)

    dδϕ_dr   = diff_r(δϕ, dr, Nr)
    d2δϕ_dr2 = diff_r(dδϕ_dr, dr, Nr)

    lhs0 = zeros(Float64, Nr)
    lhs2 = zeros(Float64, Nr)
    for i=2:Nr
        lhs0[i] = d2ϕ_dr2[i] + 2/r[i]*dϕ_dr[i] 
        lhs2[i] = d2δϕ_dr2[i] + 2/r[i]*dδϕ_dr[i] - β^2/r[i]^2*δϕ[i]
    end

    for i = Nco:Nco+10
        println(r[i], " ", T1[i], " ", T2[i], " ", ξr[i], " ", ξt[i])
    end

    plot(r, T1)

end

main()

#A = [1.0 -2.0; -1.0 1.0]
#b = [1.0, -2.0]

#    display(A)
#    display(b)
#    c = A \ b
#    display(A*c)
#    display(inv(A)*A)
