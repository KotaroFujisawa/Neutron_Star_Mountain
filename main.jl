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


    S = zeros(Float64, Nr)
    #gravitational potential
    for i=1:Nr
        for j=1:Nr
            if r[i] ≤ r[j]
               S[j] = 4π*r[j]*ρ[j]
            else
               S[j] = 4π*r[j]^2/r[i]*ρ[j]
            end
        end
        ϕ[i] =-integ_simpson(S, Nr, dr)
    end
    plot(r, ϕ, xlims=(0.0,1.0))
    K = (-ϕ[1] + ϕ[Nr]) / (Γ / (Γ - 1))

    for i=1:Nr
        p[i]   = K*ρ[i]^Γ
        cs2[i] = K*Γ*ρ[i]^(Γ-1)
    end
###end make background


###perturbations

    ft = zeros(Float64, Nr)

    for i=1:Nr
        ft[i] = 0.1*( 1 - r[i])^2
    end
    Nstep = 10
    for step = 1:Nstep

        for i=1:Nr
            δρ[i] = -(ρ[i]*δϕ[i] - ft[i] * r[i]) / cs2[i]
        end
        δρ[Nr] = 0.0

        for i=1:Nr
            δϕ[i] = 0.0
        end

    end
    
    A = [1.0 -2.0; -1.0 1.0]
    b = [1.0, -2.0]

    plot(r, δρ, xlims=(0.0,1.0))

#    display(A)
#    display(b)
#    c = A \ b
#    display(A*c)
#    display(inv(A)*A)
end

main()