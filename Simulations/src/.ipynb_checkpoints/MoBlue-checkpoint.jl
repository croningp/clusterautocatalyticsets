function make_MoBlueCRS_brokendown(k_f::Float64, k_d_stable::Float64, dimerization_ratio::Float64, mo6_enhance::Float64 = 1.0, mo36_enhance::Float64 = 1.0, mo154_enhance::Float64 = 1.0, mo132_enhance::Float64 = 1.0,volume::Float64 = 1.0, T::Float64 = 300.0, R::Float64= 8.3144598, MoMass::Int64= 95)
    """
    Arguements:
        k_f - The generic forward reaction constant
        stable_backward - The reverse reactation constant for stable molecules
        dimerization_ratio - The ratio between the reaction constant for edge dimers and corner dimers formation
        mo36_enhance - The degree to which Mo36 enhances the formation of pentamers
        wheel_enhancement_multiplier - The slope of the line which increases the formation rate of wheels as they approach completion
        ball_growth_multipler - The slope of the line which increases the formation rate of balls as they approach completion
        volume - Reaction Volume
        T - Temperature
        R - Gas constant
        MoMass - Mass of Mo
    """
    
    ########################################################################################
    # This is a GIANT FUNCTION! It could be broken down into smaller bits but it is not yet.
    # All of the "physics" of the Mo Blue is contained within this function, it specificies
    # which reactions are allowed, and how the rate constants are assigned. 
    ########################################################################################
    
    ######### Stable molecules that won't undergo the reverse reaction
    stable_molecules = ["Mo1", "Mo2", "edgeMo2", "Mo6", "Mo36","Mo36*(Mo6)14_(Mo2)28_(Mo1)14","(Mo6)14_(Mo2)28_(Mo1)14", "(Mo6)12_(edgeMo2)30"]
    
    ######### All the fragments of the various assemblies
    penta_fragments = make_penta_fragments()
    wheel_fragments = make_wheel_fragments()
    ball_fragments = make_ball_fragments()
    
    ######### Group all the different fragments and make a dictionary to map molecule strings to indicies
    all_molecules = Array{String, 1}(undef, 0)
    append!(all_molecules, penta_fragments)
    append!(all_molecules, ["edgeMo2"])
    append!(all_molecules, wheel_fragments)
    append!(all_molecules, ball_fragments)
    
    molecule_dict = Dict{String,Int64}()
    for (n, m) in enumerate(all_molecules)
        molecule_dict[m] = n
    end
    
    ######## Initialize a reaction list and rxn id number ############
    reaction_list = Array{Reaction,1}(undef,0)
    rID = 1
    
    reaction_list = generate_Mo1_Mo8_synthesis_reactions(reaction_list, all_molecules, molecule_dict, k_f, k_d_stable, dimerization_ratio, 0.0, T,R,volume)
    reaction_list = generate_Mo8_Mo36_synthesis_reactions(reaction_list, all_molecules, molecule_dict, k_f, k_d_stable, T,R,volume)
    reaction_list = generate_Mo154_synthesis_reactions(reaction_list, all_molecules, molecule_dict, k_f, k_d_stable, mo154_enhance, T, R, volume)
    reaction_list = generate_Mo132_synthesis_reactions(reaction_list, all_molecules, molecule_dict, k_f, k_d_stable, mo132_enhance, T,R, volume)
    reaction_list = generate_Mo36_templating_Mo6_synthesis_reactions(reaction_list, all_molecules, molecule_dict, k_f, k_d_stable, mo36_enhance, T,R, volume)
    reaction_list = generate_Mo6_templating_Mo6_synthesis_reactions(reaction_list, all_molecules, molecule_dict, k_f, k_d_stable, mo6_enhance, T,R,volume)
    
    reaction_list = remove_duplicate_reactants(reaction_list)
    
    parameters = Dict(
    :k_f => k_f,
    :stable_backward=> k_d_stable,
    :dimerization_ratio => dimerization_ratio,
    :mo6_enhance => mo6_enhance,
    :mo36_enhance => mo36_enhance,
    :mo154_enhance => mo154_enhance,
    :mo132_enhance => mo132_enhance,
    :volume => volume,
    :T => T,
    :R => R
    )
    
    MoBlue_CRS = CRS(all_molecules, molecule_dict, reaction_list, parameters)
   
    for rID in 1:length(MoBlue_CRS.reaction_list)
        rxn = MoBlue_CRS.reaction_list[rID]
        if length(rxn.reactants) != length(rxn.react_coef)
            println(format_reaction_str(MoBlue_CRS, rID))
            println("Reactant's dont match, RXN ID: ",rID)
            println(rxn.reactants, ", ", rxn.react_coef)
        elseif length(rxn.products) != length(rxn.prod_coef)
            println(format_reaction_str(MoBlue_CRS, rID))
            println("Products's dont match, RXN ID: ",rID)
            println(rxn.reactants, ", ", rxn.react_coef)
        end
    end
    
    return MoBlue_CRS
end
function format_reaction_str(CRS, rID)
    rxn = CRS.reaction_list[rID]
    react_IDs = rxn.reactants
    prod_IDs = rxn.products
    rx_str = ""
    for r in react_IDs
        rx_str= rx_str*CRS.molecule_list[r]* " + "
    end
    rx_str = rx_str[1:end-2]
    rx_str = rx_str*"--> "
    for p in prod_IDs
        rx_str= rx_str*CRS.molecule_list[p]* " + "
    end
    rx_str = rx_str[1:end-2]
    return rx_str
end

function count_mo_num(molecule)
    """Count the number of Mo in a molecule string """
    n = 0 # init to zero
    
    if occursin("_", molecule) # See if molecule has multiple parts or if it's a pentamer or smaller
        if occursin(")_(", molecule)
            n += 1
        end
        components = split(molecule, "_")
        for c in components
            if occursin("*", c)
                c = split(c, "*")[2]
                n += 36
            end
            size= 0
            motif, multiplier_string = split(c, ")")
            if occursin("+", motif)
                motif, add = split(motif, "+")
                size += parse(Int64, add)
            end

            if multiplier_string != ""
                multiplier = parse(Int64, multiplier_string)
            else
                multiplier =1
            end

            size += parse(Int64,split(motif, "o")[2])

            n +=  size*multiplier
        end
        
    elseif occursin("*", molecule) # See if the molecule is bound to a template or has a hanger
        addition = split(molecule, "*")[1]
        molecule = split(molecule, "*")[2]
        if occursin("hanger", addition) #Its got a hanger
            n +=1
        elseif occursin("Mo6", addition) #It's templating 
            n += 6
        elseif occursin("Mo36", addition) #It's templating 
            n += 36
        end
        
        size = split(molecule, "o")[2] # Count up the rest
        if occursin("+", size)
            
            size, adds = split(size, "+")
            n += parse(Int64,adds)
        end
        
        size = parse(Int64, size)
        n += size
    else # Otherwise it's a pentamer or smaller
        size = split(molecule, "o")[2]
        if occursin("+", size)
            
            size, adds = split(size, "+")
            n += parse(Int64,adds)
        end
        
        size = parse(Int64, size)
        n += size
    end
        
    return n 
end
##############################################################################################################
function generate_Mo132_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d_stable, mo132_enhance, T,R, volume)
    rID = length(reaction_list) + 1
    edgeMo2_id= molecule_dict["edgeMo2"]
    Mo6_id = molecule_dict["Mo6"]
    ball_frags = molecule_list[occursin.("_(edgeMo2", molecule_list)]
    
    ## Make association between pentamer and edgeMo2
    bonded_id = molecule_dict["(Mo6)1_(edgeMo2)1"]
    forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num("edgeMo2"), count_mo_num("Mo6"))
    forward_rxn = Reaction(rID, [edgeMo2_id,Mo6_id], [1,1], [bonded_id], [1],[],[], "STD", forward_k)
    push!(reaction_list, forward_rxn)
    rID += 1
    backward_k = 1.0 # intermediates are always unstable
    backward_rxn = Reaction(rID, [bonded_id], [1], [edgeMo2_id,Mo6_id], [1,1],[],[], "STD", backward_k)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    for r1 in ball_frags
        r1_id = molecule_dict[r1]
        p = add_pentamer_ball(r1)
        
        if p in ball_frags
            p_id = molecule_dict[p]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo6"], "   Product: ", prod_f)
            k_ball_f = calculate_ball_addition_rate(k_f, mo132_enhance, r1, "Mo6")
            forward_k = bimolecular_coef(k_ball_f,T,R,volume,count_mo_num(r1), count_mo_num("Mo6"))
            if p == "(Mo6)12_(edgeMo2)30"
                backward_k = k_d_stable
            else
                backward_k = 1.0
            end
            forward_rxn = Reaction(rID, [r1_id,Mo6_id], [1,1], [p_id], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1
            backward_rxn = Reaction(rID, [p_id], [1], [r1_id,Mo6_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1   
        end

        p2 = add_dimer_ball(r1)
        if p2 in ball_frags
            p2_id = molecule_dict[p2]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "edgeMo2"], "   Product: ", prod_f)
            k_ball_f = calculate_ball_addition_rate(k_f, mo132_enhance, r1, "edgeMo2")
            forward_k = bimolecular_coef(k_ball_f,T,R,volume,count_mo_num(r1), count_mo_num("edgeMo2"))
            if p == "(Mo6)12_(edgeMo2)30"
                backward_k = k_d_stable
            else
                backward_k = 1.0
            end
            forward_rxn = Reaction(rID, [r1_id,edgeMo2_id], [1,1], [p2_id], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1
            backward_rxn = Reaction(rID, [p2_id], [1], [r1_id,edgeMo2_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
        end           
    end

   return reaction_list
end
##############################################################################################################
function count_ball_components(m)
    """
    Counts the components of MoBlue Keggin "Balls" and fragments of balls
    Arguement: m - molecule string
    returns: n_pent - the number of pentamer units in the ball fragment
             n_dimer - the number of dimer units in the ball fragment
    """
    n_pent = 0
    n_dimer = 0
    pent_string = split(m,"_")[1]
    dimer_string = split(m,"_")[2]
    
    n_pent = parse(Int64,split(pent_string, ")")[2])
    n_dimer = parse(Int64,split(dimer_string, ")")[2])
    
    return n_pent, n_dimer
end

function sigmoid(x,L,k, x0) 
    return L/(1  + exp(-k*(x-x0)))
end

function calculate_ball_addition_rate(k, multiplier,m, addition)
    n_max_dimer = [5,9,12,14,18,20,23,25,27,29,30,30]
    n_min_dimer = [0,1,3,5,7,10,12,15,18,21,25,30 ]
    
    n_pent, n_dimer = count_ball_components(m)
    max_d = n_max_dimer[n_pent]
    min_d = n_min_dimer[n_pent]
    k_f = k
    if addition == "edgeMo2"
        k_f = k_f*(1+ multiplier*(max_d - n_dimer))
    elseif addition == "Mo6"
        k_f = k_f*(1+multiplier*(n_dimer- min_d))
    end
    
    return k_f
end

function build_ball_frag(n_pent, n_dimers)
    """
    Builds a string represented a ball fragment with a specificed number of pentamers and dimers
    Arguements: n_pent - number of pentamers in fragment
                n_dimers - number of dimers in fragment
    Returns: molecule string representing fragment
    """
    return "(Mo6)"*string(n_pent)*"_(edgeMo2)"*string(n_dimers) 
end

function add_dimer_ball(m)
    """
    Adds a dimer to a ball fragment, returns the updated molecule string
    """
    n_pent, n_dimer = count_ball_components(m)
    return build_ball_frag(n_pent, n_dimer+1)
end
function add_pentamer_ball(m)
    """
    Adds a pentamer to a ball fragment, returns the updated molecule string 
    """
    n_pent, n_dimer = count_ball_components(m)
    return build_ball_frag(n_pent+1, n_dimer)
end
function make_ball_fragments()
    """
    Makes all possible ball fragments allowed in the simulation
    returns an array of molecule strings representing the fragments
    """
    n_p = 1:12
    n_max_dimer = [5,9,12,14,18,20,23,25,27,29,30,30]
    n_min_dimer = [0,1,3,5,7,10,12,15,18,21,25,30 ]

    ball_fragments = Array{String,1}(undef,0)
    for p in n_p
        num_dimers = [i for i in n_min_dimer[p]:n_max_dimer[p]]
        for n in num_dimers
            if !(p == 1 && n == 0)
                push!(ball_fragments,"(Mo6)"*string(p)*"_(edgeMo2)"*string(n) )
            end

        end
    end
    return ball_fragments
end
##############################################################################################################
function generate_Mo154_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d_stable, mo154_enhance, T,R, volume)
    rID = length(reaction_list) + 1
    Mo2_id= molecule_dict["Mo2"]
    Mo1_id = molecule_dict["Mo1"]
    Mo6_id = molecule_dict["Mo6"]
    wheel_frags = molecule_list[occursin.("Mo36*(Mo6)", molecule_list)]
    #####################################################################################
    ####### Make the initial association reactions and dissociate Mo36 at the end #######
    r1 = "Mo36"
    r2 = "Mo6"
    r1_id = molecule_dict["Mo36"]
    r2_id = Mo6_id
    p1 = "Mo36*(Mo6)1_(Mo2)0_(Mo1)0"
    p1_id = molecule_dict[p1]
    k_f_wheel = calculate_wheel_stability(k_f, mo154_enhance, p1)
    forward_k = bimolecular_coef(k_f_wheel,T,R,volume,count_mo_num(r1), count_mo_num(r2))
    backward_k = 1.0
    forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p1_id], [1],[],[], "STD", forward_k)
    push!(reaction_list, forward_rxn)
    rID += 1
    backward_rxn = Reaction(rID, [p1_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
    push!(reaction_list, backward_rxn)
    rID += 1
    # Dissociate Quickly
    r1 = "Mo36*(Mo6)14_(Mo2)28_(Mo1)14"
    r1_id = molecule_dict[r1]
    p1 = "Mo36"
    p1_id = molecule_dict[p1]
    p2 = "(Mo6)14_(Mo2)28_(Mo1)14"
    p2_id = molecule_dict[p2]
    forward_k = 1000.0
    forward_rxn = Reaction(rID, [r1_id], [1], [p1_id, p2_id], [1, 1],[],[], "STD", forward_k)
    push!(reaction_list, forward_rxn)
    rID += 1
    
    # Next make all the reactions to assemble the wheel
    for r1 in wheel_frags
        r1_id = molecule_dict[r1]
        # Add a pentamer to this fragment
        p1 = add_pent_wheel(r1)
        # Check if that product is possible 
        if p1 in molecule_list 
            # Whats the product and record the reaction
            p1_id = molecule_dict[p1]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo6"], "   Product: ", prod_f)
            k_f_wheel = calculate_wheel_stability(k_f, mo154_enhance, p1)
            forward_k = bimolecular_coef(k_f_wheel,T,R,volume,count_mo_num(r1), count_mo_num("Mo6"))
            if p1 == "Mo36*(Mo6)14_(Mo2)28_(Mo1)14"
                backward_k = k_d_stable
            else
                backward_k = 1.0
            end
            
            forward_rxn = Reaction(rID, [r1_id,Mo6_id], [1,1], [p1_id], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1
            backward_rxn = Reaction(rID, [p1_id], [1], [r1_id,Mo6_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
        end
        # A a monomer to the wheel and see if its a possible product that hasn't been recorded yet
        p2 = add_monomer_wheel(r1)
        # Check if that product is possible 
        if p2 in molecule_list 
            # Whats the product and record the reaction
            p2_id = molecule_dict[p2]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo6"], "   Product: ", prod_f)
            k_f_wheel = calculate_wheel_stability(k_f, mo154_enhance, p2)
            forward_k = bimolecular_coef(k_f_wheel,T,R,volume,count_mo_num(r1), count_mo_num("Mo1"))
            if p1 == "Mo36*(Mo6)14_(Mo2)28_(Mo1)14"
                backward_k = k_d_stable
            else
                backward_k = 1.0
            end
            
            forward_rxn = Reaction(rID, [r1_id,Mo1_id], [1,1], [p2_id], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1
            backward_rxn = Reaction(rID, [p2_id], [1], [r1_id,Mo1_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
        end
        # Add a dimer to the wheel fragment and see if it's possible and whether it's been recorded 
        p3 = add_dimer_wheel(r1)
        # Check if that product is possible 
        if p3 in molecule_list 
            # Whats the product and record the reaction
            p3_id = molecule_dict[p3]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo6"], "   Product: ", prod_f)
            k_f_wheel = calculate_wheel_stability(k_f, mo154_enhance, p3)
            forward_k = bimolecular_coef(k_f_wheel,T,R,volume,count_mo_num(r1), count_mo_num("Mo2"))
            if p1 == "Mo36*(Mo6)14_(Mo2)28_(Mo1)14"
                backward_k = k_d_stable
            else
                backward_k = 1.0
            end
            
            forward_rxn = Reaction(rID, [r1_id,Mo2_id], [1,1], [p3_id], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1
            backward_rxn = Reaction(rID, [p3_id], [1], [r1_id,Mo2_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
        end
    end
    return reaction_list
end
##############################################################################################################
function build_wheel_frag(n_pent, n_dimer, n_monomer)
    """ Make a wheel molecule string with the given component counts """
    frag = "Mo36*(Mo6)"*string(n_pent)*"_(Mo2)"*string(n_dimer)*"_(Mo1)"*string(n_monomer)
    return frag
end
function calculate_wheel_stability(k, multiplier, m)
    """Calculate the wheel forward reaction rate for a given molecule m, base rate k, and slope multiplier """
    total_fragments = 14 + 28 + 14
    n_pent, n_dimer, n_monomer = count_wheel_fragments(m)
    k_f = k*(1+ multiplier*((n_pent + n_dimer + n_monomer)/total_fragments) )
    return k_f 
end
function count_wheel_fragments(m)
    """Count the number of components in the wheel """
    n_pent = 0
    n_dimer = 0
    n_monomer = 0
    pent_string = split(m,"_")[1]
    dimer_string = split(m,"_")[2]
    monomer_string = split(m,"_")[3]
    
    n_pent = parse(Int64,split(pent_string, ")")[2])
    n_dimer = parse(Int64,split(dimer_string, ")")[2])
    n_monomer = parse(Int64,split(monomer_string, ")")[2])
    
    return n_pent, n_dimer, n_monomer
end
function add_pent_wheel(m)
    """Add a pentamer to a wheel """
    pent, dimer, monomer = count_wheel_fragments(m)
    new_m = build_wheel_frag(pent + 1, dimer, monomer)
    return new_m
end
function add_dimer_wheel(m)
    """Add a corner bonded dimer to the wheel """
    pent, dimer, monomer = count_wheel_fragments(m)
    new_m = build_wheel_frag(pent, dimer+1, monomer)
    return new_m
end
function add_monomer_wheel(m)
    """Add a monomer to the wheel """
    pent, dimer, monomer = count_wheel_fragments(m)
    new_m = build_wheel_frag(pent, dimer, monomer+1)
    return new_m
end
function make_wheel_fragments()
    """Generate a list of all possible wheel fragments. This asserts that the
    wheel must be built in a particular order """
    frags = Array{String,1}(undef,0)
    n_dimer = 0
    n_monomer = 0
    for pent in 1:14
        
        push!(frags, build_wheel_frag(pent, n_dimer, n_monomer))
        n_monomer += 1
        push!(frags, build_wheel_frag(pent, n_dimer, n_monomer))
        n_dimer += 1
        push!(frags, build_wheel_frag(pent, n_dimer, n_monomer))
        n_dimer += 1
        push!(frags, build_wheel_frag(pent, n_dimer, n_monomer))

    end
    mo_assembled = split(frags[end], "*")[2]
    push!(frags, mo_assembled)
    return frags 
end
##############################################################################################################
function generate_Mo128_synthesis_reactions(reactions, molecule_list, molecule_dict, k_f, k_d, Mo128_enhance)
    
    return reactions, molecules 
end
##############################################################################################################
##############################################################################################################
function generate_Mo1_Mo8_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d_stable, dimerization_ratio, trimer_stability, T,R,volume)
    #############################################################
    ############ YOU HAVEN"T ADDED TRIMER STABILITY! ############
    #############################################################
    
    # Initialize some values
    rID = length(reaction_list) +1
    monomer_id = molecule_dict["Mo1"]
    molecule_sizes = count_mo_num.(molecule_list)
    # Get all the possible molecules between Mo1 to hanger*Mo6+3 (Mo10) exclude aggregates like (Mo6)1_(Mo2)1
    stable_molecules= ["Mo1", "Mo2", "edgeMo2", "Mo6"]
    relevant_molecules = molecule_list[molecule_sizes .<= 10]
    relevant_molecules = relevant_molecules[.!occursin.("_", relevant_molecules)]
    relevant_molecules = relevant_molecules[.!occursin.("Mo6*", relevant_molecules)]
    relevant_molecules = relevant_molecules[relevant_molecules .!= "edgeMo2"]
    ####################################################
    ##### Handle the formation of edgeMo2 manually #####
    r1 = "Mo1"
    r1_id = molecule_dict[r1]
    p1 = "edgeMo2"
    p1_id = molecule_dict[p1]
    forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r1))  
    backward_k = k_d_stable 
    forward_rxn = Reaction(rID, [r1_id], [2], [p1_id], [1],[],[], "STD", dimerization_ratio*forward_k)
    push!(reaction_list, forward_rxn)
    rID += 1
    backward_rxn = Reaction(rID, [p1_id], [1], [r1_id], [2],[],[], "STD", backward_k)
    push!(reaction_list, backward_rxn)
    rID += 1
    ###################################################################
    ##### Check all pairwise interactions between these molecules #####
    n = length(relevant_molecules)
    for i in 1:n, j in i:n
        r1 = relevant_molecules[i]
        r2 = relevant_molecules[j]
        p, extra_mo1 = merge_penta_fragments(r1,r2)
        if p in relevant_molecules
            # Get the IDs of the relevant molecules
            r1_id = molecule_dict[r1]
            r2_id = molecule_dict[r2]
            p_id = molecule_dict[p]
            # Calculate the forward and backward reaction rate constant
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r2))
            if p in stable_molecules
                backward_k = k_d_stable
            else
                backward_k = 1.0
            end
            if extra_mo1 == 0 # If the reaction doesn't have extra Mo1 produced include the reverse reaction
                if r1 != r2
                    forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id], [1],[],[], "STD", forward_k)
                    push!(reaction_list, forward_rxn)
                    rID += 1
                    backward_rxn = Reaction(rID, [p_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
                    push!(reaction_list, backward_rxn)
                    rID += 1
                elseif r1 == r2
                    forward_rxn = Reaction(rID, [r1_id], [2], [p_id], [1],[],[], "STD", forward_k)
                    push!(reaction_list, forward_rxn)
                    rID += 1
                    backward_rxn = Reaction(rID, [p_id], [1], [r1_id], [2],[],[], "STD", backward_k)
                    push!(reaction_list, backward_rxn)
                    rID += 1
                end
            else #If the reaction makes extra Mo1 you don't add the reverse reaction since it's already included in the above step 
                if r1 != r2
                    forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id, monomer_id], [1, extra_mo1],[],[], "STD", forward_k)
                    push!(reaction_list, forward_rxn)
                    rID += 1
                elseif r1 == r2
                    forward_rxn = Reaction(rID, [r1_id], [2], [p_id, monomer_id], [1, extra_mo1],[],[], "STD", forward_k)
                    push!(reaction_list, forward_rxn)
                    rID += 1
                end
            end
        end
    end
    return reaction_list 
end
##############################################################################################################
##############################################################################################################
function generate_Mo3_Mo12_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d, Mo3_template_effect)
    
    
    return reactions, molecules 
end
##############################################################################################################
##############################################################################################################
function generate_Mo8_Mo36_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d_stable, T,R,volume)
    rID = length(reaction_list) +1 
    monomer_id = molecule_dict["Mo1"]
    molecule_sizes = count_mo_num.(molecule_list)
    relevant_molecules = molecule_list[molecule_sizes .< 36]
    relevant_molecules = relevant_molecules[relevant_molecules .!= "edgeMo2"]
    ####################################################
    #### Generate Stacks from pentamer combinations #### 
    pentamers = relevant_molecules[.!occursin.("_(", relevant_molecules)]
    pentamers = pentamers[.!occursin.("Mo6*", pentamers)]
    np = length(pentamers)
    for i in 1:np, j in i:np
        r1 = pentamers[i]
        r2 = pentamers[j]
        if check_if_stackable(r1,r2)
            #### Forward and backward reactions
            p, extra_mo1 = stack(r1,r2)
            r1_id = molecule_dict[r1]
            r2_id = molecule_dict[r2]
            p_id = molecule_dict[p]
            
            # Calculate the forward and backward reaction rate constant
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r2))
            backward_k = 1.0 # Stacks are never stable 
            
            if extra_mo1 == 0 # If the reaction doesn't have extra Mo1 produced include the reverse reaction
                forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id], [1],[],[], "STD", forward_k)
                push!(reaction_list, forward_rxn)
                rID += 1
                backward_rxn = Reaction(rID, [p_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
                push!(reaction_list, backward_rxn)
                rID += 1
            else #If the reaction makes extra Mo1 you don't add the reverse reaction since it's already included in the above step 
                forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id, monomer_id], [1, extra_mo1],[],[], "STD", forward_k)
                push!(reaction_list, forward_rxn)
                rID += 1
            end
        end
    end
    ###############################################
    ####### Generate Mo36 by merging stacks #######
    stacks = relevant_molecules[occursin.(")_(", relevant_molecules)]
    # println(stacks)
    ns = length(stacks)
    for i in 1:ns
        r1 = stacks[i]
        
        #### Try adding monomer to stack ####
        possible_additions = add_monomer_to_stack_all(r1)
        if length(possible_additions)> 0
            r2 = "Mo1"
            r2_id = molecule_dict["Mo1"]
            r1_id = molecule_dict[r1]
            for p in possible_additions
                p_id = molecule_dict[p]
                forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r2))
                backward_k = 1.0 # Stacks are never stable 
                forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id], [1],[],[], "STD", forward_k)
                push!(reaction_list, forward_rxn)
                rID += 1
                backward_rxn = Reaction(rID, [p_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
                push!(reaction_list, backward_rxn)
            end
        end
    end
    for i in 1:ns, j in i:ns
        #### Check stack combinations ####
        r1 = stacks[i]
        r2 = stacks[j]
        # println(i," , ", j)
        # println(stacks[i]," , ", stacks[j])
        # println(r1," , ", r2)
        stack_possible, extra_mo1 = can_merge_stacks(r1,r2)
        if stack_possible
            #### Forward and backward reactions
            r1_id = molecule_dict[r1]
            r2_id = molecule_dict[r2]
            p = "Mo36"
            p_id = molecule_dict["Mo36"] 
            # Calculate the forward and backward reaction rate constant
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r2))
            backward_k = k_d_stable # Mo36 is always stable
            
            if extra_mo1 == 0 # If the reaction doesn't have extra Mo1 produced include the reverse reaction
                forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id], [1],[],[], "STD", forward_k)
                push!(reaction_list, forward_rxn)
                rID += 1
                backward_rxn = Reaction(rID, [p_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
                push!(reaction_list, backward_rxn)
                rID += 1
            else #If the reaction makes extra Mo1 you don't add the reverse reaction since it's already included in the above step 
                forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id, monomer_id], [1, extra_mo1],[],[], "STD", forward_k)
                push!(reaction_list, forward_rxn)
                rID += 1
            end
        end
    end
    return reaction_list
end
##############################################################################################################
function generate_Mo6_templating_Mo6_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d_stable, Mo6_enhance,T,R, volume)
    rID = length(reaction_list)+1
    Mo6_id = molecule_dict["Mo6"]
    monomer_id = molecule_dict["Mo1"]
    relevant_molecules = ["Mo1","Mo2", "Mo3", "Mo4","Mo5"]
    
    ### Generate all attachement reactions (very slow detachement)
    nm = length(relevant_molecules)
    for i in 1:nm
        r1 = relevant_molecules[i]
        r2 = "Mo6"
        p = "Mo6*"*r1
        
        r1_id = molecule_dict[r1]
        r2_id = molecule_dict[r2]
        p_id = molecule_dict[p]
        
        forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r2))
        backward_k = 1.0
        forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id], [1],[],[], "STD", Mo6_enhance*forward_k)
        push!(reaction_list, forward_rxn)
        rID += 1
        backward_rxn = Reaction(rID, [p_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
        push!(reaction_list, backward_rxn)
        rID += 1
        ### Generate all forward reactions (basically polymerization), no backward reactions
        p2 = "Mo6*"*add_mo_penta(r1)
        p2_id = molecule_dict[p2]
        forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(p), count_mo_num("Mo1"))
        forward_rxn = Reaction(rID, [p_id,monomer_id], [1,1], [p2_id], [1],[],[], "STD", Mo6_enhance*forward_k)
        push!(reaction_list, forward_rxn)
        rID += 1
    end
    ### Dissociate Mo6*Mo6 quickly
    r1 = "Mo6*Mo6"
    r1_id = molecule_dict[r1]
    p1 = "Mo6"
    p1_id = molecule_dict[p1]
    backward_k = 1.0
    backward_rxn = Reaction(rID, [r1_id], [1], [p1_id], [2],[],[], "STD", 1000*backward_k)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    return reaction_list 
end
##############################################################################################################
function generate_Mo36_templating_Mo6_synthesis_reactions(reaction_list, molecule_list, molecule_dict, k_f, k_d_stable, Mo36_enhance, T,R,volume)
    rID = length(reaction_list)+1
    Mo36_id = molecule_dict["Mo36"]
    monomer_id = molecule_dict["Mo1"]
    
    relevant_molecules = ["Mo1","Mo2", "Mo3", "Mo4","Mo5"]
    ### Generate all attachement reactions (very slow detachement)
    nm = length(relevant_molecules)
    for i in 1:nm
        r1 = relevant_molecules[i]
        r2 = "Mo36"
        p = "Mo36*"*r1
        
        r1_id = molecule_dict[r1]
        r2_id = molecule_dict[r2]
        p_id = molecule_dict[p]
      
        forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(r1), count_mo_num(r2))
        backward_k = 1.0 
        forward_rxn = Reaction(rID, [r1_id,r2_id], [1,1], [p_id], [1],[],[], "STD", Mo36_enhance*forward_k)
        push!(reaction_list, forward_rxn)
        rID += 1
        backward_rxn = Reaction(rID, [p_id], [1], [r1_id,r2_id], [1,1],[],[], "STD", backward_k)
        push!(reaction_list, backward_rxn)
        rID += 1
        ### Generate all forward reactions (basically polymerization), no backward reactions
        p2 = "Mo36*"*add_mo_penta(r1)
        p2_id = molecule_dict[p2]
        forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(p), count_mo_num("Mo1"))
        forward_rxn = Reaction(rID, [p_id,monomer_id], [1,1], [p2_id], [1],[],[], "STD", Mo36_enhance*forward_k)
        push!(reaction_list, forward_rxn)
        rID += 1
    end
    ### Dissociate Mo36*Mo6 quickly
    r1 = "Mo36*Mo6"
    r1_id = molecule_dict[r1]
    p1 = "Mo36"
    p1_id = molecule_dict[p1]
    p2 = "Mo6"
    p2_id = molecule_dict[p2]
    backward_rxn = Reaction(rID, [r1_id], [1], [p1_id, p2_id], [1,1],[],[], "STD", 1000*1.0)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    return reaction_list 
end
##############################################################################################################
function remove_duplicate_reactants(reaction_list)
    nr = length(reaction_list)
    for i in 1:nr
        rxn = reaction_list[i]
        reactants = rxn.reactants
        if length(reactants)> 1
            if reactants[1] == reactants[2]
                new_reaction = Reaction(rxn.index, [reactants[1]], [2], rxn.products, rxn.prod_coef, rxn.catalysts, rxn.cat_coef, rxn.propensity, rxn.rate_constant)
                reaction_list[i] = new_reaction
            end
        end
    end
    return reaction_list
end
##############################################################################################################
function build_penta_fragment(core,additions,hanger)
    """Build a pentamer fragment with the listed number of mo """
    m ="NULL"
    if hanger>=2
        return m
    elseif (core < 6 && hanger > 0)
        return m
    elseif additions >3
        return m
    elseif (core <3 && additions >=1)
        return m
    elseif (core <4 && additions >=2)
        return m
    elseif (core <6 && additions >=3)
        return m
    else
        m = "Mo"*string(core)
        if additions > 0
            m= m*"+"*string(additions)
        end
        if hanger != 0
            m = "hanger*"*m
        end
        return m
    end
end
function count_components_penta(m)
    """Count how many Mo are in a pentamer fragment """
    mo = 0
    additions = 0
    hanger = 0
    if m == "NULL"
        hanger = 0
        mo = 0
        additions = 0
    else
        front = ""
        if occursin("*", m)
            front = split(m, "*")[1]
            m = split(m, "*")[2]
        end
        hanger = Int(occursin("hanger",front))
        if occursin("+", m)
            addition_string = split(m, "+")[2]
            additions = parse(Int64,addition_string)
        else
            additions = 0
        end
        mo_string = split(m, "+")[1]
        mo = parse(Int64, split(mo_string, "o")[2])
    end
    return mo,additions, hanger
end
function add_hanger_penta(m)
    """ Add a hanger to a pentamer fragment """
    mo, additions, hanger = count_components_penta(m)
    return build_penta_fragment(mo,additions, hanger+1)
end
function add_addition_penta(m)
    """ Add an addition the pentamer fragment """
    mo, additions, hanger = count_components_penta(m)
    return build_penta_fragment(mo,additions+1, hanger)
end
function add_mo_penta(m)
    """ Add a Mo to the core of a pentamer fragment """
    mo, additions, hanger = count_components_penta(m)
    return build_penta_fragment(mo+1,additions, hanger)
end

function merge_penta_fragments(m1,m2)
    """ Merge two pentamer fragments and see if there are extra monomers produced """
    ##################################################################################
    # Pentamer fragments can merge if either have a core of less than 6 (we assume the 
    # stability of the Mo6 pentamer prevents it), but an Mo6 might pick up hangers and
    # extras from the fragments
    ##################################################################################
    extra = 0 # extra monomers
    m_new = "NULL"
    mo1,add1,hanger1 = count_components_penta(m1)
    mo2,add2,hanger2 = count_components_penta(m2)
    mo_sum = sum([mo1, add1, hanger1, mo2, add2, hanger2])
    if mo1 ==6 && mo2== 6
        m_new = "NULL"
        extra = 0
    elseif mo_sum >= 10
        extra = mo_sum -10 
        m_new = build_penta_fragment(6, 3, 1)
    else 
        extra = 0
        additions = 0
        hangers = 0
        core = minimum([mo_sum, 6])
        if core > 5
            additions = minimum([mo_sum-6,3])
        else
            additions = 0
            hangers = 0
        end
        if additions == 3 && mo_sum == 10
            hangers = 1
        else
            hangers = 0
        end
            
        m_new = build_penta_fragment(core, additions, hangers)
    end
    if any([occursin("Mo36*", m1), occursin("Mo36*", m2)])
        m_new = "NULL"
        extra = 0
    end
    if m_new == "hanger*Mo6+3" && any([m1,m2] .== "hanger*Mo6+3")
        m_new = "NULL"
        extra = 0
    end
    return m_new, extra
end
function make_penta_fragments()
    """ Generate all possible pentamer fragments based on rules """
    frags = Array{String,1}(undef, 0)
    for i in 1:6
        ## Add catalyst assemblies 
        cataylzed_frag = "Mo36*"*build_penta_fragment(i,0,0)
        push!(frags,cataylzed_frag)
        cataylzed_frag = "Mo6*"*build_penta_fragment(i,0,0)
        push!(frags,cataylzed_frag)
        for j in 0:3
            m1 = build_penta_fragment(i,j, 0)
            if m1 != "NULL"
                push!(frags, m1)
            end
            m2 = build_penta_fragment(i,j,1)
            if m2 != "NULL"
                push!(frags,m2)
            end
        end
    end
    # push!(frags, "edgeMo2")
    stacks = Array{String,1}(undef,0)
    for p in frags
        for q in frags
            if check_if_stackable(p,q)
                push!(stacks,stack(p,q)[1])
            end
        end
    end
    append!(frags, stacks)
    push!(frags, "Mo36")
    return frags
end

function check_if_stackable(m1,m2)
    """ See if two pentamer fragments can be stacked """
    stack = false
    m1_hanger = occursin("hanger",m1)
    m2_hanger = occursin("hanger",m2)
    m1_mo36 = occursin("Mo36", m1)
    m2_mo36 = occursin("Mo36", m2)
    if any([m1_hanger, m2_hanger]) && !any([m1_mo36, m2_mo36])
        if (occursin("Mo6", m1) && occursin("Mo6", m2))
            stack = true
        end
    end
    return stack 
end

function stack(m1,m2)
    """ Stack two pentamer fragments """
    extra = 0
    m1_hanger = occursin("hanger", m1)
    m2_hanger = occursin("hanger", m2)
    if m1_hanger
        m1 = split(m1, "*")[2]
    end
    if m2_hanger
        m2 = split(m2, "*")[2]
    end
    stack = "("*m1*")_("*m2*")"
    if m1_hanger && m2_hanger
        extra = 1
    end
    return stack,extra
end

function can_add_monomer_to_stack(m)
    """See if a monomer can be added to a stack in any way """
    possible =false
    left,right  = split(m, ")_(")
    if occursin("+", left)
        left = split(left, "+")[2]
        left_adds = parse(Int64,split(left, ")")[1])
    else
        left_adds= 0 
    end
    if left_adds < 3
        possible = true
    end
    if occursin("+", right)
        right = split(right, "+")[2]
        right_adds = parse(Int64,split(right, ")")[1])
    else
        right_adds= 0 
    end
    if right_adds <3
        possible = true
    end
    return possible
end

function can_remove_monomer_from_stack(m)
    """ See if a monomer can come off a stack""" 
    possible =false
    left,right  = split(m, ")_(")
    if occursin("+", left)
        left = split(left, "+")[2]
        left_adds = parse(Int64,split(left, ")")[1])
    else
        left_adds= 0 
    end
    if left_adds > 0
        possible = true
    end
    if occursin("+", right)
        right = split(right, "+")[2]
        right_adds = parse(Int64,split(right, ")")[1])
    else
        right_adds= 0 
    end
    if right_adds > 0
        possible = true
    end
    return possible
end

function add_monomer_to_stack_all(m)
    """ Add monomer to stack, check all possible outcomes """
    all_possible = Array{String,1}(undef,0)
    left,right  = split(m, ")_(")
    left = split(left, "(")[2]
    right = split(right, ")")[1]
    core_left,add_left,hanger = count_components_penta(left)
    core_right,add_right,hanger = count_components_penta(right)
    new_left = build_penta_fragment(core_left, add_left + 1, 1)
    if new_left != "NULL"
        push!(all_possible,stack(new_left, right)[1])
    end
    new_right = build_penta_fragment(core_right, add_right+1,1)
    if new_right !="NULL"
        push!(all_possible,stack(left,new_right)[1])
    end
    return all_possible
end

function can_merge_stacks(m1,m2)
   """See if two stacks can be merged and if extra monomers are produced in the process"""
    extra = 0 
    possible = false
    left1, right1 = split(m1,")_(")
    left2, right2 = split(m2,")_(")
    if occursin("+",left1)
        left1 = split(left1, "+")[2]
        left1_adds = parse(Int64,split(left1, ")")[1])
    else
        left1_adds = 0
    end
    if occursin("+",right1)
        right1 = split(right1, "+")[2]
        right1_adds = parse(Int64,split(right1, ")")[1])
    else
        right1_adds = 0
    end
    if occursin("+",left2)
        left2_adds = parse(Int64,split(left2, "+")[2])
    else
        left2_adds = 0
    end
    if occursin("+",right2)
        right2 = split(right2, "+")[2]
        right2_adds = parse(Int64,split(right2, ")")[1])
    else
        right2_adds = 0
    end
    if sum([left1_adds, left2_adds, right1_adds, right2_adds]) >= 10
        if all([left1_adds, left2_adds, right1_adds, right2_adds] .> 1)
            possible = true
            extra = sum([left1_adds, left2_adds, right1_adds, right2_adds]) - 10
        end
    end
    return possible, extra 
end


