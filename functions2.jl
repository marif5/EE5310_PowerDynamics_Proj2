using JuMP, Ipopt

# ===================== Power blance fucntion for Gen Injections ===========================#

function Powerbalance(Bus_V, Bus_ang, Y_adm, Pd, Qd, wind_MW)
                                                                #= Test Values =#
                                                                # Bus_V      = results_ACPF.Bus_V
                                                                # Bus_ang    = results_ACPF.Bus_ang
                                                                # i = 0;
                                                                # wind_MW = zeros(73,1)
                                                                # Y_adm = Y
                                                                # Pd = Bus_Pd
                                                                # Qd = Bus_Qd
                                                                # for w in length(wind_data[:,1])
                                                                #     i = i + 1;
                                                                #     wind_MW[w,1] = wind_data[i,5];
                                                                # end
    Vc      = Bus_V.*exp.(im*Bus_ang)
    Ic      = Y_adm*Vc
    Sc      = Vc.*conj.(Ic)
    P_b     = real(Sc) 
    Q_b     = imag(Sc)

    Pload  = Pd;
    Qload  = Qd;

    Pgen_inj = zeros(73,1);
    Qgen_inj = zeros(73,1);

    Pgen_inj = -wind_MW./100 + P_b + Pload
    Qgen_inj =  Q_b + Qload
    
    return Pgen_inj, Qgen_inj
end

# ===================== Classical Generator initialization fucntion ===========================#

function initialize_cgs(p, q, v, delta_bus, xdp,bno)
                                                            #= Test Values =#
                                                            # p          = Pg_inj[classical_gen];
                                                            # q         = Qg_inj[classical_gen];
                                                            # v          = V_b[classical_gen];
                                                            # delta_bus     = ang_b[classical_gen];
                                                            # xdp       = classical_data[:,5];
                                                            # bno       = length(classical_data[:,1])
    model = Model(Ipopt.Optimizer)
    @variable(model, (delta_bus[j]-pi/4) <= delta_gen[j=1:bno] <= (delta_bus[j]+pi/4))
    @variable(model, 0.5 <= Ep[j = 1:bno] <= 1.5)
    # constraints
    for ind in 1:bno
        @constraint(model, p[ind] == (v[ind]*Ep[ind]/xdp[ind])*sin(delta_gen[ind]-delta_bus[ind]))
        @constraint(model, q[ind] == (v[ind]*Ep[ind]/xdp[ind])*cos(delta_gen[ind]-delta_bus[ind]) - (v[ind]^2)/xdp[ind])
    end
    @objective(model, Min, 1)
    optimize!(model)

    # make sure the solution is valid
    println(termination_status(model))
    println(Int(termination_status(model)))
    if (Int(termination_status(model)) != 4)
        @assert 1 == 2
    end

    return value.(Ep), value.(delta_gen)
end

# ===================== Classical Generator initialization fucntion ===========================#
# ## Q3 -----------------------------------------------------------------------
function unreduced_jacobian(V,ang,Y)
    # Calculating J11, J12, J21, J22
    J11 = zeros(Buses, Buses)
    J12 = zeros(Float64, Buses, Buses)
    J21 = zeros(Float64, Buses, Buses)
    J22 = zeros(Float64, Buses, Buses)
    J = zeros(146,146)
    sorted_PV_dict = keys(sort(PV_dict))
    slack_index = keys(sort(slack_dict)) 
        # Calculate diagonal elements
        for k in 1:Buses, n in 1:Buses
            if k != n
                # diagonal
                J11[k, k] += abs(V[k]) * abs(V[n]) * abs(Y[k, n]) * sin(angle(Y[k, n]) - ang[k] + ang[n])
                J12[k, k] = 2 * abs(V[k]) * abs(Y[k, k]) * cos(angle(Y[k, k])) + (sum(abs(V[n]) * abs(Y[k, n]) * cos(angle(Y[k, n]) - ang[k] + ang[n]) for n in 1:Buses if k != n))
                J21[k, k] += abs(V[k]) * abs(V[n]) * abs(Y[k, n]) * cos(angle(Y[k, n]) - ang[k] + ang[n])
                J22[k, k] = -2 * abs(V[k]) * abs(Y[k, k]) * sin(angle(Y[k, k])) - (sum(abs(V[n]) * abs(Y[k, n]) * sin(angle(Y[k, n]) - ang[k] + ang[n]) for n in 1:Buses if k != n))
                        
                # off-diagonal
                J11[k, n] = -(abs(V[k]) * abs(V[n]) * abs(Y[k, n]) * sin(angle(Y[k, n]) - ang[k] + ang[n]))
                J12[k, n] = abs(V[k]) * abs(Y[k, n]) * cos(angle(Y[k, n]) - ang[k] + ang[n])
                J21[k, n] = -abs(V[k]) * abs(V[n]) * abs(Y[k, n]) * cos(angle(Y[k, n]) - ang[k] + ang[n])
                J22[k, n] = -abs(V[k]) * abs(Y[k, n]) * sin(angle(Y[k, n]) - ang[k] + ang[n])
             end
    
        end
        J = [J11 J12; J21 J22]
        return J
    end

function generator_jacobian(id, Ep, delta_g, V, delta_b, xdp)
    Ep = Ep[1:end];
    V = V[1:end];
    delta_b = delta_b[1:end];
    delta_g = delta_g[1:end];
    cg_buses = id[1:end];
    no_cgens = length(cg_buses) #number of classical generators
    
    Jac_gen = zeros(146,146);
    J11_entries = zeros(no_cgens);
    J11 = zeros(73,73);
    J12_entries = zeros(no_cgens);
    J12 = zeros(73,73);
    J21_entries = zeros(no_cgens);
    J21 = zeros(73,73);
    J22_entries = zeros(no_cgens);
    J22 = zeros(73,73);

        for ii in 1:no_cgens
            J11_entries[ii] = (Ep[ii]*sin(delta_g[ii] - delta_b[ii]))/xdp[ii];
            J12_entries[ii] = -V[ii]*Ep[ii]*cos(delta_g[ii] - delta_b[ii])/xdp[ii];
            J21_entries[ii] = (Ep[ii]*cos(delta_g[ii] - delta_b[ii])/xdp[ii]) - 2*V[ii]/xdp[ii];            
            J22_entries[ii] = V[ii]*Ep[ii]*sin(delta_g[ii] - delta_b[ii])/xdp[ii];
            
            J11[cg_buses[ii],cg_buses[ii]] = J11_entries[ii];
            J12[cg_buses[ii],cg_buses[ii]] = J12_entries[ii];
            J21[cg_buses[ii],cg_buses[ii]] = J21_entries[ii];
            J22[cg_buses[ii],cg_buses[ii]] = J22_entries[ii];
        end
    
    Jac_gen = [J21 J22; J11 J12];
end

function compute_PQcg(Ep,V,xdp,delta_g,delta_b)
    Pcg = zeros(73,1)
    Qcg = zeros(73,1)
    i = 0
    for j in classical_gen
        i = i + 1
        Pcg[j] = (Ep[i]*V[i]/xdp[i])*sin(delta_g[i] - delta_b[i])
        Qcg[j] = (Ep[i]*V[i]/xdp[i])*cos(delta_g[i] - delta_b[i])-V[i]^2/xdp[i]
    end
end

function recompute_residual(Bus_V, Bus_ang, Y_adm, Pd, Qd, wind_MW)
    Vc      = Bus_V.*exp.(im*Bus_ang)
    Ic      = Y_adm*Vc
    Sc      = Vc.*conj.(Ic)
    P_b     = real(Sc) 
    Q_b     = imag(Sc)

    Pload  = Pd;
    Qload  = Qd;

    Pgen_inj = zeros(73,1);
    Qgen_inj = zeros(73,1);
    
                Ep = df_cg.Ep
                xdp = df_cg_para.xdp
                delta_g = df_cg.delta_g
                delta_b = df_cg.delta_b
                V = Bus_V

                Pcg = zeros(73,1)
                Qcg = zeros(73,1)
                i = 0
                for j in classical_gen
                    i = i + 1
                    Pcg[j] = (Ep[i]*V[i]/xdp[i])*sin(delta_g[i] - delta_b[i]) 
                    Qcg[j] = (Ep[i]*V[i]/xdp[i])*cos(delta_g[i] - delta_b[i])-V[i]^2/xdp[i] 
                end

    Pgen_inj = -wind_MW./100 + P_b + Pload -Pcg
    Qgen_inj =  Q_b + Qload - Qcg
                
    return Pgen_inj, Qgen_inj
    
end

function newton(int_algbpf_Jac, int_residual, V, delta, Y, Bus_Pd, Bus_Qd, wind_MW, classical_gen,df_cg)
   
    al_pf_jac  = int_algbpf_Jac
    residual   = int_residual
    
    # while norm(residual,2) > 1e-5
    for f in 1:2
            Vd_new          = [V; delta] - al_pf_jac\residual
            V               = Vd_new[1:73]              #V
            delta           = Vd_new[74:end]            #deltas
            Pinj, Qinj      = recompute_residual(V, ang, Y, Bus_Pd, Bus_Qd, wind_MW)
            
            Vcg             = V[classical_gen]
            delta_cg        = delta[classical_gen]            
            Unred_Jac       = unreduced_jacobian(V, delta, Y)
            G_Jac           = generator_jacobian(df_cg.Bus_id, df_cg.Ep, df_cg.delta_g, Vcg, delta_cg, df_cg_para.xdp);
 
            al_pf_jac        = Unred_Jac + G_Jac;
            residual        = [Pinj; Qinj];
            println(norm(residual),2) 
          
    end
    return V, delta
        
end

function fwd_euler(X, df_cg, df_cg_para, Pg_inj, V_nr , delta_nr, classical_gen)
    #stating coefficents:
       
        H      = df_cg_para.H
        Pm     = Pg_inj[classical_gen]
        D      = df_cg_para.D
        xdp    = df_cg_para.xdp
        Vb     = V_nr[classical_gen]
        delta_b= delta_nr[classical_gen]
        omega_0 = 2*pi*60*ones(length(df_cg.omega_g),1)
        Ep_0        = df_cg.Ep
               
# Forward Euler formation  
delta_g = zeros(19,1) 
omega_g = zeros(19,1)  
delta_g_dot = zeros(19,1) 
omega_g_dot = zeros(19,1)
    h = 0.02    
    times = 0.02:h:0.04
    for i in 1:length(times)
        #stating states:
        for gen in 1:19
            delta_g[gen,1]     = X[gen,1]  
            omega_g[gen,1]     = X[gen,2]  
            delta_g_dot[gen,1] = omega_g[gen,1]
            omega_g_dot[gen,1] = (omega_0[gen]/2*H[gen])*(Pm[gen] - D[gen]*omega_g[gen] - (Vb[gen]*Ep[gen]/xdp[gen])*sin(delta_g[gen]- delta_b[gen]))
        end
      
        #state equations
        # delta_g_dot = omega_g
        # omega_g_dot = (omega_0/2*H[i])*(Pm[i] - D[i]*omega_g[i] - (Vb[i]*Ep[i]/xdp[i])*sin(delta_g- delta_b[i]))

       fX = [delta_g_dot omega_g_dot]
        X  = X + h* fX
    end    
 
    return times, X      
end
