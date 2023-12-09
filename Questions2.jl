# importing libraries

using SymPy, LinearAlgebra, DataFrames, CSV, Plots, SymPy

# Including ACPF, fucntions and Gen data files

include("ACPF.jl")
include("GenData.jl")
include("functions2.jl")

## ========================== Q1: ACPF and Active and reactive power Gen injections ================================
# %%

# A - Solving ACPF for initial Bus Voltage and delta

ang, V, conv,iter = ACPF()
results_ACPF = DataFrame(
    Bus_ang = ang, 
    Bus_V = V, 
    iteration = iter
    )
    x = 1:length(results_ACPF.Bus_V)
    y = results_ACPF.Bus_V
    plot(x, y, label="ACPF Bus Voltages", xlabel="Buses", ylabel="Voltage Magnitude (V)", title="ACPF Results")
    savefig("Q1.png")
    x = 1:length(results_ACPF.Bus_ang)
    y = results_ACPF.Bus_ang
    plot(x, y, label="ACPF Phase Angles", xlabel="Buses", ylabel="Phase Angles", title="ACPF Results")
    savefig("Q11.png")

    # B - Finding initial Genertor injections using power balance

P_inj   = zeros(73,1)
Q_inj   = zeros(73,1)
P_wind = zeros(73,1);

V_b      = results_ACPF.Bus_V;
ang_b    = results_ACPF.Bus_ang;
    i = 0;
    for w in length(wind_data[:,1])
        i = i + 1;
        P_wind[w,1] = wind_data[i,5];
    end

Pg_inj, Qg_inj = Powerbalance(V_b, ang_b, Y, Bus_Pd, Bus_Qd, P_wind)
    x = 1:length(Pg_inj)
    y = Pg_inj
    plot(x, y, label="ACPF Generator Power Injections", xlabel="Buses", ylabel="Active Power (pu)", title="ACPF Results")
    savefig("Q12.png")
    x = 1:length(Qg_inj)
    y = Qg_inj
    plot(x, y, label="ACPF Generator Power Injections", xlabel="Buses", ylabel="Reactive Power (pu)", title="ACPF Results")
    savefig("Q13.png")

#ClassicalGen Data
#{bus ID, w0, H, D, xd'}
# Q2
classical_gen = Int.(classical_data[:,1]);
p_cg          = Pg_inj[classical_gen];
q_cg          = Qg_inj[classical_gen];
v_cg          = V_b[classical_gen];
delta_cg      = ang_b[classical_gen];
xdp_cg        = classical_data[:,5];
bno_cg        = length(classical_data[:,1]);


## ========================== Q2: Initialize the Classical Generators ================================
# %%

Ep,delta_gen = initialize_cgs(p_cg, q_cg, v_cg, delta_cg, xdp_cg,bno_cg)
x = 1:length(Ep)
y = Ep
plot(x, y, label="EMF voltage", xlabel="Buses", ylabel="E'", title="Classical Generator Initializations")
savefig("Q21.png")
x = 1:length(delta_gen)
y = delta_gen
plot(x, y, label="Generator angle", xlabel="Buses", ylabel="Generator angle", title="Classical Generator Initializations")
savefig("Q22.png")


# Q3    
    df_cg = DataFrame( #ACPF results and initialized variables
    Bus_id = classical_gen,
    V = results_ACPF.Bus_V[classical_gen],    
    delta_b = results_ACPF.Bus_ang[classical_gen],
    Pg = p_cg,
    Qg = q_cg,
    Ep = Ep,
    delta_g = delta_gen,
    omega_g = zeros(Float64,19)
    )

    df_cg_para = DataFrame( #data frame to keep classical generator data
    w0 = classical_data[:,2], #w0
    H = classical_data[:,3],
    D = classical_data[:,4],
    xdp = classical_data[:,5]
    )
## ========================== Q3: Simulate the system (Forward Euler) ============================
# %%

# A - Algebraic Power Flow Jacobian 
  
    Unreduced_Jac = unreduced_jacobian(results_ACPF.Bus_V, results_ACPF.Bus_ang, Y)
    Gen_Jac       = generator_jacobian(df_cg.Bus_id, df_cg.Ep, df_cg.delta_g, df_cg.V, df_cg.delta_b,df_cg_para.xdp);
    int_algb_pf_jac   =  Unreduced_Jac+Gen_Jac 
    int_residual      = [Pg_inj ; Qg_inj]
    
    V_nr, delta_nr = newton(int_algb_pf_jac, int_residual, results_ACPF.Bus_V , results_ACPF.Bus_ang , Y, Bus_Pd , Bus_Qd, P_wind, classical_gen, df_cg)
# %%
# x = 1:length(int_algb_pf_jac)
# y = int_algb_pf_jac
# spy(y,  title="Generator Jacobian")
# savefig("Q32.png")



    X = [df_cg.delta_g zeros(19,1)]
    times, state_out = fwd_euler(X, df_cg, df_cg_para, Pg_inj, V_nr , delta_nr, classical_gen)
    plot(state_out[:,1])
    plot!(state_out[:,2])


    x = 1:length(state_out[:,1])
    y = state_out[:,1]
    z = state_out[:,2]
    plot(x, y, label="Generator Angle", title="Forward Euler results")
    plot!(x, z, label="Generator Angular frequency", title="Forward Euler results")
    savefig("Q33.png")
    