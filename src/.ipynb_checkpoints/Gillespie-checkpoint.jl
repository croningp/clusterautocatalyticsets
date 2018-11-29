 #module MultilevelGenesis

#############################################
###### Dependencies #########################
#############################################
using Random
#using Plots
using DataFrames
using JSON
using CSV
using Combinatorics

include("Data.jl")
include("Model_CRS.jl")
include("MoBlue.jl")
#############################################
###### Rection Functions ####################
#############################################

function CalculatePropensities(concentrations::Array{Float64,1}, CRS::CRS) 
    ########## NOT UNIT TESTED! ##########
    n_rxn = length(CRS.reaction_list)::Int64
    PropensityVec = zeros(n_rxn)
    rID = 0
    reaction_list = CRS.reaction_list::Array{Reaction,1}
    for rxn in CRS.reaction_list
        Ap = 0.0
        if rxn.propensity == "STD"
            Ap = STD_propensity(rxn, concentrations)::Float64
        end
        rID += 1
        #println("reaction_list ID: ", rID, "    Propensity: ", Ap)
        PropensityVec[rID] = Ap
    end
    
    return PropensityVec
end

function PickReactionID(PropensityVec::Array{Float64,1})
    Ap_tot = sum(PropensityVec)
    dice_roll = Ap_tot*rand()
    checkpoint = 0.0
    
    rID = 0
    #println("Dice Roll: ", dice_roll)
    while checkpoint < dice_roll
        rID += 1
        checkpoint += PropensityVec[rID]
        #println("rID: ", rID, "     Checkpoint: ", checkpoint)
    end
    #println("Returning: ", rID, "   Rxn Propensity: ", PropensityVec[rID])
   return rID 
end

function ExecuteReaction(concentrations::Array{Float64,1}, CRS::CRS, rID::Int64)
    ########## NOT UNIT TESTED! ##########
    rxn = CRS.reaction_list[rID]
    for (i, r) in enumerate(rxn.reactants)
        concentrations[r] -= rxn.react_coef[i]
    end
    for (j, p) in enumerate(rxn.products)
        concentrations[p] += rxn.prod_coef[j]
    end
    if any(x->x<0, concentrations)
        println(CRS.reaction_list[rID])
        for (i, r) in enumerate(rxn.reactants)
            concentrations[r] += rxn.react_coef[i]
        end
        for (j, p) in enumerate(rxn.products)
            concentrations[p] -= rxn.prod_coef[j]
        end
        println("Concentration of Reactants: ",concentrations[rxn.reactants])
        println("Rxn Propensity based on concentrations: ", STD_propensity(rxn, concentrations))
        error("Concentrations went negative on reaction $rID")
        
    end
    
    return concentrations
end

function STD_propensity(rxn::Reaction, concentrations::Array{Float64,1})
    ########## NOT UNIT TESTED! ##########
    Ap  = rxn.rate_constant
    for (ri,r) in enumerate(rxn.reactants)
        rconc = concentrations[r]
        r_order = rxn.react_coef[ri]
        Ap = Ap*(prod([(rconc-i) for i in 0:(r_order-1)])) 
    end
    cat_enhance = 0.0
    for (ci,c) in enumerate(rxn.catalysts)
        cconc = concentrations[c]
        cat_enhance += cconc*rxn.cat_coef[ci]
    end
    Ap = Ap*(1 + cat_enhance)
    return Ap
end

# function sort_reactions_by_propensity(concentrations, myCRS)
#     println(myCRS.reaction_list[1:15])
#     propensities = CalculatePropensities(concentrations, myCRS)
#     println(propensities)
#     sorted_indices = sortperm(propensities, rev = true)
#     println(sorted_indices)
#     myCRS.reaction_list[:] = myCRS.reaction_list[sorted_indices]
#     println(myCRS.reaction_list[1:15])
#     return myCRS
# end

#############################################
##### Stochastic Evolution Functions ########
#############################################

 
function Run_Binary_Ligation(nA::Int64, nB::Int64, max_L::Int64, kf::Float64, kb::Float64, volume::Float64, tau_max::Float64, tau_freq::Float64 = 0.1, T::Float64 = 300.0, R::Float64= 8.3144598, mA::Int64 = 100, mB::Int64 = 100)
    ### Unique Experiment number required
    
    ## Check for the correct output files and directories
    param_file = "../data/Binary_Ligation_Run_Parameters.csv"
    save_dir = "../data/Binary_Ligation_timeseries"
    
    if !isdir(save_dir)
        mkdir(save_dir)
    end
    
        
    ## Initialize Time values
    tau = 0.0
    freq_count = 0.0
    
    binary_CRS = generate_binary_ligation_CRS(max_L, kf, kb, volume, T)
    ### Redo output initialization using DataFrames
    output_DF = InitializeOutput(binary_CRS)
               
    nmolecules = length(binary_CRS.molecule_list)
    concentrations= zeros(nmolecules)
    
    Aindex = binary_CRS.molecule_dict["A"]
    Bindex = binary_CRS.molecule_dict["B"]
    concentrations[Aindex] = nA
    concentrations[Bindex] = nB
        
    #Initialize Propensities
    propensities = CalculatePropensities(concentrations, binary_CRS)

    ####### Main LOOP #######
    while tau < tau_max
        # Pick Reaction
        rID = PickReactionID(propensities)
        # Execute Reaction
        concentrations = ExecuteReaction(concentrations, binary_CRS, rID)
        # Calculate Propensities
        propensities = CalculatePropensities(concentrations, binary_CRS)
        # Record Data
        if tau >= freq_count
            output_DF[Symbol(freq_count)] = concentrations[:] 
            freq_count += tau_freq
            println(tau)
        end

        #Update Time
        Ap_tot = sum(propensities)
        tau -= (log(rand())/Ap_tot) 
    end

    # Save time series and params
    save_name, parameter_df = generate_output_data(binary_CRS, save_dir)#
    println(parameter_df)
    if !isfile(param_file)
        CSV.write(param_file, parameter_df)
    else
        runs_df = CSV.read(param_file)
        runs_df = [runs_df; parameter_df]
    end
    
    CSV.write(save_name, output_DF)
    return concentrations
end

function Run_MoBlue(Mo1_mass::Int64, total_time::Float64, k_f::Float64, stable_backward::Float64, dimerization_ratio::Float64, mo36_enhance::Float64, wheel_enhancement_multiplier::Float64= 1.0,ball_growth_multipler::Float64 = 1.0,out_count::Float64 = 100.0, volume::Float64 = 1.0, T::Float64 = 300.0, R::Float64= 8.3144598)::Array{Float64,1}
    
    ### Unique Experiment number required
    tau_freq = total_time/out_count
    ## Check for the correct output files and directories
    param_file = "../data/MoBlue_Run_Parameters.csv"
    save_dir = "../data/MoBlue_timeseries"
    
    if !isdir(save_dir)
        mkdir(save_dir)
    end
    
        
    ## Initialize Time values
    tau = 0.0
    freq_count = 0.0
    
    MoBlueCRS = make_MoBlueCRS(k_f, stable_backward, dimerization_ratio, mo36_enhance,wheel_enhancement_multiplier,ball_growth_multipler, volume, T, R)::CRS
    ### Redo output initialization using DataFrames
    output_DF = InitializeOutput(MoBlueCRS)
               
    nmolecules = length(MoBlueCRS.molecule_list)::Int64
    concentrations= zeros(nmolecules)::Array{Float64,1}
    
    Monomer_index = MoBlueCRS.molecule_dict["Mo1"]
    concentrations[Monomer_index] = Mo1_mass
    
    #Initialize Propensities
    propensities = CalculatePropensities(concentrations, MoBlueCRS)

    ####### Main LOOP #######
    while tau < total_time
        # Pick Reaction
        rID = PickReactionID(propensities)::Int64
        #println("rID is : ", rID, "     PropensityVec: ", propensities[rID] )
        #println("Picked reaction is: ", MoBlueCRS.reaction_list[rID])
        #println("Reactants: ", MoBlueCRS.molecule_list[MoBlueCRS.reaction_list[rID].reactants])
        # println("Products: ", MoBlueCRS.molecule_list[MoBlueCRS.reaction_list[rID].products])
        #println("Concentrations: ", concentrations[MoBlueCRS.reaction_list[rID].reactants])
        #println("Propensity: ", STD_propensity(MoBlueCRS.reaction_list[rID],concentrations))
        
        # Execute Reaction
        concentrations = ExecuteReaction(concentrations, MoBlueCRS, rID)::Array{Float64,1}
        
        # Calculate Propensities
        propensities = CalculatePropensities(concentrations, MoBlueCRS)::Array{Float64,1}
        # Record Data
        if tau >= freq_count
            # Calculate mass and check conservation...
            current_mass = calculate_mass(concentrations, MoBlueCRS)::Float64
            if current_mass != Mo1_mass
                println(MoBlueCRS.reaction_list[rID])
                println(MoBlueCRS.molecule_dict[MoBlueCRS.reaction_list[rID].reactants],MoBlueCRS.molecule_list[MoBlueCRS.reaction_list[rID].products] )
                error("Conservation of Mass not maintained.. Exiting Simulation")
            end
            output_DF[Symbol(freq_count)] = concentrations[:] 
            freq_count += tau_freq
            #MoBlueCRS = sort_reactions_by_propensity(concentrations, MoBlueCRS)
            
            
            sort!(MoBlueCRS.reaction_list, by=x -> STD_propensity(x, concentrations), rev = true)
            propensities = CalculatePropensities(concentrations, MoBlueCRS)::Array{Float64,1}
            #println(tau)
        end

        #Update Time
        Ap_tot = sum(propensities)
        tau -= (log(rand())/Ap_tot) 
    end

    # Save time series and params
    save_name, parameter_df = generate_output_data(MoBlueCRS, save_dir)#
    #println(parameter_df)
    if !isfile(param_file)
        #println("Not Finding File")
        CSV.write(param_file, parameter_df)
    else
        CSV.write(param_file, parameter_df; append = true)
    end
    CSV.write(save_name, output_DF)
    return concentrations
end

function calculate_mass(concentrations, CRS)
    mass = 0
    nmolecules= length(CRS.molecule_list)
    for i in 1:nmolecules
        M = count_mo_num(CRS.molecule_list[i])
        mass += M*concentrations[i]
    end
    return mass
end



#############################################
#end #### module #############################
#############################################