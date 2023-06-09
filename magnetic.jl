using DifferentialEquations
using Dierckx
using Plots
using PyCall
@pyimport matplotlib.pyplot as plt

function make_mf()

    Bp0 = 1.0e12
    Bp(r) = Bp0*(1.0 - 0.99*(r/r_max)^4)
    β = 1*(1+2)
    function source(du, u, param, r)
        #u[1] -> m(r), u[2] -> p(r)
        Br = u[1]
        du[1] = -2*Br / r + β*Bp(r)/r
    end
    param = [0.0]
    r_min = 1.0e3
    r_max = 1.0e6
    rspan = (r_min, r_max)
    Br0 = 1.0e12
    u0 = [Br0]
    prob = DifferentialEquations.ODEProblem(source, u0, rspan, param)
    bg_mf = solve(prob, Tsit5(), reltol = 1.0e-12, abstol = 1.0e-12)

    nr = 100
    r_g = range(r_min, r_max, nr)

    fig, ax = plt.subplots()

    ax.plot(bg_mf.t, bg_mf.u)
    ax.grid(true)
    ax.set_ylabel("Br")

    Bp_g = zeros(Float64, nr)
    for i=1:nr
        Bp_g[i] = Bp(r_g[i])
    end
    ax.plot(r_g, Bp_g)

    ax.set_yscale("log")
    
    ax.plot(r_g, Bp_g)
 
    plt.savefig("mf.pdf")



#    display(bg_mf.u)

end

make_mf()