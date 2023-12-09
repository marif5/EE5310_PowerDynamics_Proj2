# importing libraries
# %%
using SymPy, LinearAlgebra, DataFrames, CSV, Plots, SymPy

#importing MAT file data as CSVs to DataFrames

path_1 = "/Users/beyzanuraydin/Desktop/UVM/Fall2023/Energy Sys Analys/project2/Sam Project 1/BranchData.csv"
path_2 = "/Users/beyzanuraydin/Desktop/UVM/Fall2023/Energy Sys Analys/project2/Sam Project 1/BusData.csv"
path_3 = "/Users/beyzanuraydin/Desktop/UVM/Fall2023/Energy Sys Analys/project2/Sam Project 1/loadData.csv"
path_4 = "/Users/beyzanuraydin/Desktop/UVM/Fall2023/Energy Sys Analys/project2/Sam Project 1/GenData.csv"	
path_5 = "/Users/beyzanuraydin/Desktop/UVM/Fall2023/Energy Sys Analys/project2/Sam Project 1/wind_Data.csv"	

BranchData = CSV.File(path_1) |> DataFrame
BusData = CSV.File(path_2) |> DataFrame
loadData = CSV.File(path_3) |> DataFrame
GenData = CSV.File(path_4) |> DataFrame	
WindData = CSV.File(path_5) |> DataFrame		
Lines = size(BranchData,1)	
 

BranchData.Z = BranchData.R + im*BranchData.X
BranchData.Y = 1 ./BranchData.Z
BranchData.ont .= parse.(Complex{Float64}, BranchData.ont)

loadData.Pd = loadData.Pd./100
theta = acos.(loadData.Pf)	
loadData.Qd = tan.(theta).*loadData.Pd	


# ======== #

output_dict = Dict{Int, Float64}()
for i in 1:length(GenData.bus)
    key = GenData.bus[i]
    value = GenData.Vreg[i]
    if !haskey(output_dict, key)
        output_dict[key] = value
    end
end	

#For Load Pd 	
Dic_1 = Dict{Int, Float64}()
for i in 1:length(loadData.bus)
    key = loadData.bus[i]
    value = loadData.Pd[i]
	 Dic_1[key] = value
end

# for Load Qd 
Dic_2 = Dict{Int, Float64}()
for i in 1:length(loadData.bus)
    key = loadData.bus[i]
    value = loadData.Qd[i]
	 Dic_2[key] = value
end

# For Generator Pg 
Dic_3  = Dict{Int, Float64}()
	for i in 1:length(GenData.Pg)
		key = GenData.bus[i]
		value = GenData.Pg[i]
		if haskey(Dic_3, key)
			Dic_3[key] += value
		else
			Dic_3[key] = value
		end
	end	
		
# For Generator Qg,max and min 
Dic_4 = Dict{Int, Tuple{Float64, Float64}}()

for i in 1:length(GenData.bus)
    key = GenData.bus[i]
    value_1 = GenData.Qmax[i]
    value_2 = GenData.Qmin[i]
    if haskey(Dic_4, key)
        Dic_4[key] = (Dic_4[key][1] + value_1, Dic_4[key][2] + value_2)
    else
        Dic_4[key] = (value_1, value_2)
    end
end

#For wind generation
	
Bus = 1:73	
Bus = collect(Bus)		
Bus_Vmag = [get(output_dict, bus, 0) == 0 ? 1 : get(output_dict, bus, 0) for bus in Bus]
Bus_Pd  = [get(Dic_1, bus, 0) == 0 ? 0 : get(Dic_1, bus, 0) for bus in Bus]
Bus_Qd  = [get(Dic_2, bus, 0) == 0 ? 0 : get(Dic_2, bus, 0) for bus in Bus]
Bus_Pg =  [get(Dic_3, bus, 0) == 0 ? 0 : get(Dic_3, bus, 0) for bus in Bus]
Bus_qmax = [get(Dic_4, bus, (0.0, 0.0))[1] for bus in Bus]	
Bus_qmin = 	[get(Dic_4, bus, (0.0, 0.0))[2] for bus in Bus]
Bus_Qg = zeros(Float64, 73)
Bus_delta  = zeros(Float64, 73)
Bus_type = BusData.Bus_type	
GenData_2 = DataFrame(Bus = Bus, Bus_type= Bus_type, Bus_mag = Bus_Vmag, Bus_delta = Bus_delta, Bus_Pd = Bus_Pd, Bus_Qd = Bus_Qd, Bus_Pg = Bus_Pg, Bus_Qg = Bus_Qg, Bus_qmax = Bus_qmax, Bus_qmin = Bus_qmin)	
GenData_2.Bus_Pg[16] = GenData_2.Bus_Pg[16] + 231.89
GenData_2.Bus_Pg[17] = GenData_2.Bus_Pg[17] + 262.41
GenData_2.Bus_Pg[18] = GenData_2.Bus_Pg[18] + 368.73
GenData_2.Bus_Pg[19] = GenData_2.Bus_Pg[19] + 144.11
GenData_2.Bus_Pg[20] = GenData_2.Bus_Pg[20] + 103.56
GenData_2.Bus_Pg[21] = GenData_2.Bus_Pg[21] + 195.11
GenData_2.Bus_Pg[53] = GenData_2.Bus_Pg[53] + 195.11
GenData_2.Bus_Pg[56] = GenData_2.Bus_Pg[56] + 52.33	
GenData_2.Bus_Pg = GenData_2.Bus_Pg./100	

GenData
GenData_2

####===
Buses = size(BusData,1)
    Y = zeros(Complex, Buses, Buses)
    # Non-diagonal elements
    for l in 1:Lines
        Y[BranchData[l,:fbus], BranchData[l,:tbus]] += (-BranchData[l,:Y])./(conj(BranchData[l,:ont]))
        Y[BranchData[l,:tbus], BranchData[l,:fbus]] += -BranchData[l,:Y]./(BranchData[l,:ont])
    end
    # Diagonal elements
    for b in 1:Buses, l in 1:Lines
        if (b==BranchData[l,:fbus])
			Y[b,b] += (BranchData[l,:Y] + 0.5*im*BranchData[l,:b] + BusData[b,:Bsh])./(abs.(BranchData[l,:ont]))^2
			(b==BranchData[l,:tbus])
		end
		if (b==BranchData[l,:tbus])
            Y[b,b] += (BranchData[l,:Y] + 0.5*im*BranchData[l,:b] + BusData[b,:Bsh])
        end
    end
	Y

    norm(Y)
    # ======

    # Known PV and Slack Variables 

	PV_dict = Dict()
    slack_dict = Dict()
	for (index, value_b) in enumerate(GenData_2.Bus_type)
        if value_b == 2     # PV 
            PV_dict[GenData_2.Bus[index]] = GenData_2.Bus_mag[index] 
        end
        if value_b == 3     # Slack 
            slack_dict[GenData_2.Bus[index]] = GenData_2.Bus_mag[index] 
        end
	end

# unknown V and theta

unknown_v_dict = Dict()
unknown_ang_dict = Dict()

for (index, value) in enumerate(GenData_2.Bus_type)
    if value == 1
        unknown_v_dict[index] = 1.0
        unknown_ang_dict[index] = 0.0
    end
end

# %%
####==================Function definitions ====================#####

####==================AC Power Flow =====-------===============#####

    function delta_y(V,ang)	
        P = zeros(Float64,Buses)
        del_P = zeros(Float64,Buses)	
        Q = zeros(Float64,Buses)
        del_Q = zeros(Float64,Buses)
        
        for i in 1:Buses
            P[i] = sum(abs(V[i]) * abs(V[j]) * abs(Y[i, j]) * cos(angle(Y[i, j]) - ang[i] + ang[j]) for j in 1:Buses) 
        end
        
        P_g_new = copy(GenData_2.Bus_Pg)
        P_d_new = copy(GenData_2.Bus_Pd)	
        del_P = P_g_new .- P_d_new .- P	
        index_to_remove = keys(slack_dict)
        P_delta = del_P[1:end .!= index_to_remove]
        
        for i in 1:Buses, j in 1:Buses
            Q[i] = -sum(abs(V[i]) * abs(V[j]) * abs(Y[i, j]) * sin(angle(Y[i, j]) - ang[i] + ang[j]) for j in 1:Buses) 
        end
        
        Q_g_new = copy(GenData_2.Bus_Qg)
        Q_g_d = copy(GenData_2.Bus_Qd)	
        del_Q = Q_g_new .- Q_g_d .- Q	
            
        sorted_v_dict = sort(unknown_v_dict)
        indices_to_keep = [parse(Int, split(string(variable), "Bus")[end]) for (variable, value) in sorted_v_dict]	
        Q_delta = del_Q[indices_to_keep]
        delta_y = vcat(P_delta, Q_delta)	
        return delta_y
        end 	

#=======#
       
function Jacob(V, ang)
# Calculating J11, J12, J21, J22
J11 = zeros(Buses, Buses)
J12 = zeros(Float64, Buses, Buses)
J21 = zeros(Float64, Buses, Buses)
J22 = zeros(Float64, Buses, Buses)
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
            J11 = J11[setdiff(1:end,keys(slack_dict)), setdiff(1:end,keys(slack_dict))]
            J_12_remove = merge!(PV_dict, slack_dict)	
            J12 = J12[1:end .∉ [slack_index], 1:end .∉ [keys(J_12_remove)]]	
            J21 = J21[1:end .∉ [keys(J_12_remove)], 1:end .∉ [slack_index]]	
            J22 = J22[1:end .∉ [keys(J_12_remove)], 1:end .∉ [keys(J_12_remove)]]	
            d =0;   # switch for decoupled powerflow

            J = [J11 d*J12; d*J21 J22]
            return J
end

function ACPF()	
    iter = 0	
    V = copy(GenData_2.Bus_mag)	
    ang = copy(GenData_2.Bus_delta)
    conv = 10
    slack_theta = keys(slack_dict)
    unknown_vari = sort(unknown_v_dict)		
    y_zero = 0
    Pinj = 0
    Qinj = 0
    while conv>1e-3 iter < 10
        y_zero  = delta_y(V,ang) 
        
        Jcalc = Jacob(V, ang)
        x_delta = inv(Jcalc)*y_zero	
        ang_delta = x_delta[1:72]
        v_delta = x_delta[73:end]
        insert!(ang_delta, first(slack_theta), 0.0)	
        ang = ang + ang_delta
        for (key, value) in zip(keys(unknown_vari), v_delta)
            V[key] = V[key] + value
    end	
        # for (key, value) in zip(keys(unknown_p), Px)
        #     Px[key] = Px[key] + value
        # end
        # for (key, value) in zip(keys(unknown_q), Qx)
        #     Qx[key] = Qx[key] + value
        # end		
            conv = norm(y_zero)
            iter = iter + 1

        
    end
        return ang, V, conv,iter
end



# # rename!(results,[:Bus_ph_ang,:Bus_Volt,:err])

# output_dec_ACPF = ACPF()
# vector1, vector2, scalar1, scalar2 = output_ACPF
# results_dec_ACPF = DataFrame(col1 = vector1, col2 = vector2, col3 = scalar1)
# rename!(results,[:Bus_ph_ang,:Bus_Volt,:err])


# Help/inspiration take from Hamza's code and GPT 3.5.
#=====================================================#