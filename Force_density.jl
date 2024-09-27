module Force_density

    function set_force_density(type_force, ρ_r, dρ_dr_r, r_c, r_o, R, A, j_0, j_1)

        # force A f_i = -Aρ ∇_i (r^2 Y_lm) Gittins et al. (2021)
        if(type_force==1)
            if(A == 0.0)
#                A = 7.9652112337329695e6
#                A = 7.777790150771358e6
                A = 7.811899165291698e6
            end
            println("Type_force==$type_force, A = $A")
            fr1(r) =-2A*r*ρ_r(r)
            ft1(r) =-A*r*ρ_r(r)
    
            return A, fr1, fr1, fr1, ft1, ft1, ft1
        end

       # force A f_i = -Aρ ∇_i (r^-3 Y_lm) Gittins et al. (2021)
        if(type_force == 2)
            if(A == 0.0)
                A = 3.6009236487675906e35
            end
            println("Type_force==$type_force, A = $A")
            fr2_co(r) = 0.0
            fr2_cr(r) = -3A*ρ_r(r) / r^4
            fr2_oc(r) = -3A*ρ_r(r) / r^4

            ft2_co(r) = 0.0
            ft2_cr(r) = A*ρ_r(r)/r^4
            ft2_oc(r) = A*ρ_r(r)/r^4

            return A, fr2_co, fr2_cr, fr2_oc, ft2_co, ft2_cr, ft2_oc
        end

        # termanl mountain Gittins et al. (2021)
        if(type_force == 3)
            if(A == 0.0)
                A = 294.4581220856396
            end
            println("Type_force==$type_force, A = $A")
            fr3(r) = -A * (2*r * ρ_r(r) + r^2 * dρ_dr_r(r))
            ft3(r) = -A * ρ_r(r) * r
            return A, fr3, fr3, fr3, ft3, ft3, ft3
        end

        # force B f_i = Bρr ∇_i Y_lm in Morales & Horowitz (2022)
        if(type_force == 4)
            if(A==0.0)
                A = 1.4086951056556168e10
            end
            println("Type_force==$type_force, A = $A")
            fr4(r) = 0.0
            ft4(r) = A*ρ_r(r)
            return A, fr4, fr4, fr4, ft4, ft4, ft4
        end

        # Morales & Horowitz (2002)
        if(type_force == 5)
            #force B f_i = Bρr ∇_i Y_lm
            if(A==0.0)
                A = 3.694291921478615e9
            end
            println("Type_force==$type_force, A = $A")
            C =  -A / (r_o-r_c)^2
            rc = 9.8e6
            function fr5(r)
                if(r >= rc)
                    return 0.0
                else
                    return 0.0
                end
            end
            function ft5(r)
                if(r >= rc)
                    return A*ρ_r(r) + C * ρ_r(r) * (r-rc)^2
                else
                    return A*ρ_r(r)    
                end
            end
            return A, fr5, fr5, fr5, ft5, ft5, ft5
        end
    
        #solenoidal force
        if(type_force == 6)
            if(A == 0.0)
                A = 0.002674122401146829
            end
            println("Type_force==$type_force, A = $A")
            fr6(r) = - 6A * r^2 * ρ_r(r) 
            ft6(r) = - 4A * r^2 * ρ_r(r)

            return A, fr6, fr6, fr6, ft6, ft6, ft6
        end

        #irrotational + solenoidal force
        if(type_force == 7)
            if(A == 0.0)
                A = 10.2674122401146829
                A = 0.0052235694170678386   # A_sol = 0.5,  A_irr = 0.5
            end
        #        A = 0.0052235694170678386   # A_sol = 0.5,  A_irr = 0.5
        #        A = 12.027155420618188      # A_sol = 0.0,  A_irr = 1.0
        #        A = 0.026163232718820248    # A_sol = 0.1,  A_irr = 0.9
        println("Type_force==$type_force, A = $A")
            A_irr = j_0
            A_sol = j_1
            fr7(r) = -A * ρ_r(r) * ( A_irr * (3 * r^2) + A_sol * (6 * r^2))
            ft7(r) = -A * ρ_r(r) * ( A_irr * (r^2)     + A_sol * (4 * r^2))

            return A, fr7, fr7, fr7, ft7, ft7, ft7
        end

        # Purely poloidal field
        if(type_force == 8)
            if(A == 0.0)
                A = 2.091119226466299e26
            end
            println("Type_force==$type_force, A = $A")
            fr8(r) = A/r^3*(6R^3*sin(π*r/R) - 6R^2*π*r*cos(π*r/R)
                        -3R*π^2*r^2*sin(π*r/R) + 3*π^3*r^3*cos(π*r/R)+2π^3*r^3
                    )*sin(π*r/R)
            ft8(r) = A/r^3*(-3R*(2R^2 - π^2*r^2)*sin(π*r/R) + π*r*(6*R^2*cos(π*r/R)+π^2*r^2)
                    )*sin(π*r/R)

            return A, fr8, fr8, fr8, ft8, ft8, ft8
        end

        # Purely poloidal field in non-barotropic 
        if(type_force == 9)
            if(A == 0.0)
                A = 4.0735347607741e24
            end
            println("Type_force==$type_force, A = $A")
            fr9(r) = A/(16π)*(1575*(r/R)^7 - 4515*(r/R)^5 + 4165(r/R)^3 - 1225(r/R))
            ft9(r) = A/(32π)*(525*(r/R)^7 - 1995*(r/R)^5 + 2695*(r/R)^3 - 1225(r/R))
            return A, fr9, fr9, fr9, ft9, ft9, ft9
        end


        # Purely poloidal field with surface current at core-crust boundary
        if(type_force == 19)
            if(A == 0.0)
                A = 2.2952271024565947e24
            end
#            F0 = 1.0
#            rho_c = 1.0
#            a(r) = 4π/3*F0*rho_c /π^5 * (
#                (3π^2*R^3*r - 6*R^5/r)*sin(π*r/R)
#                +6*π*R^4*cos(π*r/R)+π^3*R^2*r^2)
#            da_dr(r) = 4F0*rho_c/(3*π^4)*(
#                -6*R^3*π^2*sin(π*r/R) + 2*R^2*π^3*r 
#                + (6*R^5/r^2 + 3*R^3*π^2)*sin(π*r/R) 
#                + π*(-6*R^5/r + 3*R^3*π^2*r)*cos(π*r/R)/R)
#            da_dr2(r) = 4F0*rho_c/(3*π^4)*(
#                -12*R^5*sin(π*r/R)/r^3 - 6*R^2*π^3*cos(π*r/R) 
#                + 2*R^2*π^3 + 2*π*(6*R^5/r^2 + 3*R^3*π^2)*cos(π*r/R)/R 
#                - π^2*(-6*R^5/r + 3*R^3*π^2*r)*sin(π*r/R)/R^2)
            
            println("Type_force==$type_force, A = $A")
 #           println("j_0 = $j_0, j_1 = $j_1")

            fr19_co(r) = A*8*(3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) - 2*π^2*r^2*(3*R*sin(π*r/R) - π*r) - 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R) + 6*r^3*(j_0 + j_1))*sin(π*r/R)/(9*r^3)
            fr19_cr(r) = A*8*(6*R^3*sin(π*r/R) - 6*R^2*π*r*cos(π*r/R) - 3*R*π^2*r^2*sin(π*r/R) - 3*j_0*r_c^3 + 6*j_1*r^3 + 3*π^3*r^3*cos(π*r/R) + 2*π^3*r^3)*sin(π*r/R)/(9*r^3)
            fr19_oc(r) = A*8*(6*R^3*sin(π*r/R) - 6*R^2*π*r*cos(π*r/R) - 3*R*π^2*r^2*sin(π*r/R) - 3*j_0*r_c^3 - 3*j_1*r_o^3 + 3*π^3*r^3*cos(π*r/R) + 2*π^3*r^3)*sin(π*r/R)/(9*r^3)
            ft19_co(r) = A*8*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + π*r*(6*R^2*cos(π*r/R) + π^2*r^2) + 3*r^3*(j_0 + j_1))*sin(π*r/R)/(9*r^3)
            ft19_cr(r) = A*8*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r^3 + π*r*(6*R^2*cos(π*r/R) + π^2*r^2))*sin(π*r/R)/(9*r^3)
            ft19_oc(r) = A*8*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r_o^3 + π*r*(6*R^2*cos(π*r/R) + π^2*r^2))*sin(π*r/R)/(9*r^3)

            return A, fr19_co, fr19_cr, fr19_oc, ft19_co, ft19_cr, ft19_oc
        end
    
    # with current layer
        if(type_force == 20)
            if(A == 0.0)
                A = 2.2952271024565947e24
            end
        
            println("Type_force==$type_force, A = $A")

            j_2 = -1.0e-12
            h = 4.0e4

            fr10_co(r) = A*8*(3*R*(2*R^2 + pi^2*r^2)*sin(pi*r/R) + 8*h*j_2*r^3 - 2*pi^2*r^2*(3*R*sin(pi*r/R) - pi*r) - 3*pi*r*(2*R^2 - pi^2*r^2)*cos(pi*r/R))*sin(pi*r/R)/(9*r^3)
            function fr10_cr(r)
                if(r ≤ r_c+2h)
                    return (A*2*(h^2*(60*R*(2*R^2 + pi^2*r^2)*sin(pi*r/R) + pi^2*r^2*(-120*R*sin(pi*r/R) + 40*pi*r) - 60*pi*r*(2*R^2 - pi^2*r^2)*cos(pi*r/R)) + 60*j_2*r^3*(4*h^3 - 3*h*r^2 + 6*h*r*r_c - 3*h*r_c^2 + r^3 - 3*r^2*r_c + 3*r*r_c^2 - r_c^3) - j_2*(80*h^3*r^3 - 36*h*r^5 + 90*h*r^4*r_c - 60*h*r^3*r_c^2 + 6*h*r_c^5 + 10*r^6 - 36*r^5*r_c + 45*r^4*r_c^2 - 20*r^3*r_c^3 + r_c^6))*(4*R*h^2*j_2*r - 4*R*j_2*r*(h - r + r_c)^2 + h^2*pi^3*sin(pi*r/R))/(45*h^4*pi^3*r^3))
                else
                    return A*8*(15*R*(2*R^2 + pi^2*r^2)*sin(pi*r/R) - 4*h*j_2*(8*h^3 + 18*h^2*r_c + 15*h*r_c^2 + 5*r_c^3) - 10*pi^2*r^2*(3*R*sin(pi*r/R) - pi*r) - 15*pi*r*(2*R^2 - pi^2*r^2)*cos(pi*r/R))*sin(pi*r/R)/(45*r^3)
                end
            end
            fr10_oc(r) = A*8*(15*R*(2*R^2 + pi^2*r^2)*sin(pi*r/R) - 4*h*j_2*(8*h^3 + 18*h^2*r_c + 15*h*r_c^2 + 5*r_c^3) - 10*pi^2*r^2*(3*R*sin(pi*r/R) - pi*r) - 15*pi*r*(2*R^2 - pi^2*r^2)*cos(pi*r/R))*sin(pi*r/R)/(45*r^3)

            ft10_co(r) = A*8*(-3*R*(2*R^2 - pi^2*r^2)*sin(pi*r/R) + 4*h*j_2*r^3 + pi*r*(6*R^2*cos(pi*r/R) + pi^2*r^2))*sin(pi*r/R)/(9*r^3)
            function ft10_cr(r)
               if(r ≤ r_c+2h)
                    return (A*2*(20*h^2*(-3*R*(2*R^2 - pi^2*r^2)*sin(pi*r/R) + pi*r*(6*R^2*cos(pi*r/R) + pi^2*r^2)) + j_2*(80*h^3*r^3 - 36*h*r^5 + 90*h*r^4*r_c - 60*h*r^3*r_c^2 + 6*h*r_c^5 + 10*r^6 - 36*r^5*r_c + 45*r^4*r_c^2 - 20*r^3*r_c^3 + r_c^6))*(4*R*h^2*j_2*r - 4*R*j_2*r*(h - r + r_c)^2 + h^2*pi^3*sin(pi*r/R))/(45*h^4*pi^3*r^3))
                else
                    return A*8*(-15*R*(2*R^2 - pi^2*r^2)*sin(pi*r/R) + 4*h*j_2*(8*h^3 + 18*h^2*r_c + 15*h*r_c^2 + 5*r_c^3) + 5*pi*r*(6*R^2*cos(pi*r/R) + pi^2*r^2))*sin(pi*r/R)/(45*r^3)
                end
            end
            ft10_oc(r) = A*8*(-15*R*(2*R^2 - pi^2*r^2)*sin(pi*r/R) + 4*h*j_2*(8*h^3 + 18*h^2*r_c + 15*h*r_c^2 + 5*r_c^3) + 5*pi*r*(6*R^2*cos(pi*r/R) + pi^2*r^2))*sin(pi*r/R)/(45*r^3)
            return A, fr10_co, fr10_cr, fr10_oc, ft10_co, ft10_cr, ft10_oc
        end

    end

    function set_magnetic_stress(A, r_c, r_o, R, j_0, j_1)
#        println("j_0 = $j_0, j_1 = $j_1")

        m_rr_co(r) = A*4*R*(4*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + π*r*(6*R^2*cos(π*r/R) + π^2*r^2) + 3*r^3*(j_0 + j_1))^2 + (3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) - 2*π^2*r^2*(3*R*sin(π*r/R) - π*r) - 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R) + 6*r^3*(j_0 + j_1))^2)/(27*π^4*r^6)

        m_rr_cr(r) = A*4*R*(4*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r^3 + π*r*(6*R^2*cos(π*r/R) + π^2*r^2))^2 + (-3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 - 6*j_1*r^3 - 2*π^2*r^2*(-3*R*sin(π*r/R) + π*r) + 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R))^2)/(27*π^4*r^6)

        m_rr_oc(r) = A*4*R*(4*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r_o^3 + π*r*(6*R^2*cos(π*r/R) + π^2*r^2))^2 + (-3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r_o^3 - 2*π^2*r^2*(-3*R*sin(π*r/R) + π*r) + 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R))^2)/(27*π^4*r^6)

        m_rth_co(r) = A*8*R*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + π*r*(6*R^2*cos(π*r/R) + π^2*r^2)
         + 3*r^3*(j_0 + j_1))*(3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) - 2*π^2*r^2*(3*R*sin(π*r/R) - π*r) 
         - 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R) + 6*r^3*(j_0 + j_1))/(27*π^4*r^6)# / r

        m_rth_cr(r) = A*8*R*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r^3 
        + π*r*(6*R^2*cos(π*r/R) + π^2*r^2))*(3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) - 3*j_0*r_c^3 
        + 6*j_1*r^3 - 2*π^2*r^2*(3*R*sin(π*r/R) - π*r) - 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R))/(27*π^4*r^6) #/ r

        m_rth_oc(r) = A*8*R*(-3*R*(2*R^2 - π^2*r^2)*sin(π*r/R) + 3*j_0*r_c^3 + 3*j_1*r_o^3 
        + π*r*(6*R^2*cos(π*r/R) + π^2*r^2))*(3*R*(2*R^2 + π^2*r^2)*sin(π*r/R) - 3*j_0*r_c^3 
        - 3*j_1*r_o^3 - 2*π^2*r^2*(3*R*sin(π*r/R) - π*r) - 3*π*r*(2*R^2 - π^2*r^2)*cos(π*r/R))/(27*π^4*r^6)
        return m_rr_co, m_rr_cr, m_rr_oc, m_rth_co, m_rth_cr, m_rth_oc
    end
end
