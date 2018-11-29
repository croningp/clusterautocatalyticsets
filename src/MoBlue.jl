function make_MoBlueCRS(k_f::Float64, stable_backward::Float64, dimerization_ratio::Float64, mo36_enhance::Float64 = 1.0, wheel_enhancement_multiplier::Float64 = 1.0, ball_growth_multipler::Float64 = 1.0,volume::Float64 = 1.0, T::Float64 = 300.0, R::Float64= 8.3144598, MoMass::Int64= 95)
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
    
    # Stable molecules that won't undergo the reverse reaction
    stable_molecules = ["Mo1", "Mo2", "edgeMo2", "Mo6", "Mo36","Mo36*(Mo6)14_(Mo2)28_(Mo1)14","(Mo6)14_(Mo2)28_(Mo1)14", "(Mo6)12_(edgeMo2)30"]
    #stable_backward = 1E-4 
    # All the fragments of the various assemblies
    penta_fragments = make_penta_fragments()
    wheel_fragments = make_wheel_fragments()
    ball_fragments = make_ball_fragments()
    
    # Group all the different fragments and make a dictionary to map molecule strings to indicies
    all_molecules = Array{String, 1}(undef, 0)
    append!(all_molecules, penta_fragments)
    append!(all_molecules, ["edgeMo2"])
    append!(all_molecules, wheel_fragments)
    append!(all_molecules, ball_fragments)
    molecule_dict = Dict{String,Int64}()
    for (n, m) in enumerate(all_molecules)
        molecule_dict[m] = n
    end
    
    # Make a set to record all the reactions which have already been made
    reactants_set = Set()
    
    # Keep track of particular stable products
    monomer_id = molecule_dict["Mo1"]
    edge_dimer_id = molecule_dict["edgeMo2"]
    dimer_id = molecule_dict["Mo2"]
    penta_id = molecule_dict["Mo6"]
    mo36_id = molecule_dict["Mo36"]
    
    # Initialize a reaction list
    reaction_list = Array{Reaction,1}(undef,0)
    # Initialize a reaction id number
    rID = 1
    
    # Make edgeDimer formation reaction
    forward_k = bimolecular_coef(dimerization_ratio*k_f,T,R,volume,count_mo_num("Mo1"), count_mo_num("Mo1"))
    forward_rxn = Reaction(rID, [monomer_id], [2], [edge_dimer_id], [1],[],[], "STD", forward_k)
    
    #### Add the reaction to the list and increment the rID counter
    push!(reaction_list, forward_rxn)
    rID += 1
    
    #### Build the reverse reaction based on stability and add it to the list
    backward_rxn = Reaction(rID, [edge_dimer_id], [1], [monomer_id], [2],[],[], "STD", stable_backward)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    # First iterate over all the possible pentamer fragment interactions and make reactions between them
    for f1 in penta_fragments # f1 is a molecule string
         
        #### First calculate all the interactions with Monomers, because they need to be handled differently. ####
        
        m1 = molecule_dict[f1] # m1 is the index associated with f1
        # Calculate the reaction constant between f1 and a monomer
        forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num("Mo1"))
        
        if (!occursin("_", f1) && !occursin("*", f1)) # If the molecule is not a stack and doesn't have a hanger (or template)
            f1_hanger = add_hanger_penta(f1) #Add a hanger to the molecule
            if f1_hanger in all_molecules # Check if that generated molecule is valid
                p1 = molecule_dict[f1_hanger] # Get the index of the molecule
                #### Generate the forward reaction
                forward_rxn = Reaction(rID, [m1,monomer_id], [1,1], [p1], [1],[],[], "STD", forward_k)
                #### Add the reaction to the list and increment the rID counter
                push!(reaction_list, forward_rxn)
                rID += 1
                #### If the molecule is stable the reverse reaction will be VERY Slow (maybe 0?)
                if f1_hanger in stable_molecules
                    backward_k = stable_backward
                else
                    backward_k = 1.0
                end
                #### Build the reverse reaction based on stability and add it to the list
                backward_rxn = Reaction(rID, [p1], [1], [m1,monomer_id], [1,1],[],[], "STD", backward_k)
                push!(reaction_list, backward_rxn)
                rID += 1
            end
            
            f1_add = add_addition_penta(f1) # Try adding a monomer addition to the pentamer and see if it's viable
            if f1_add in all_molecules
                p1 = molecule_dict[f1_add] # p1 is the index of that molecule
                forward_rxn = Reaction(rID, [m1,monomer_id], [1,1], [p1], [1],[],[], "STD", forward_k)
                #### Add the reaction to the list and increment the rID counter
                push!(reaction_list, forward_rxn)
                rID += 1
                #### If the molecule is stable the reverse reaction will be VERY Slow (maybe 0?)
                if f1_add in stable_molecules
                    backward_k = stable_backward
                else
                    backward_k = 1.0
                end
                #### Build the reverse reaction based on stability and add it to the list
                backward_rxn = Reaction(rID, [p1], [1], [m1,monomer_id], [1,1],[],[], "STD", backward_k)
                push!(reaction_list, backward_rxn)
                rID += 1
            end
            
            f1_grow = add_mo_penta(f1) # try growing the pentamer and see if it's viable
            if f1_grow in all_molecules
                p1 = molecule_dict[f1_grow]
                #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo1"], "   Product: ", f1_grow)
                if f1 != "Mo1"
                    forward_rxn = Reaction(rID, [m1,monomer_id], [1,1], [p1], [1],[],[], "STD", forward_k)
                else
                    forward_rxn = Reaction(rID, [monomer_id], [2], [p1], [1],[],[], "STD", forward_k)
                    #println(forward_rxn)
                end
                #### Add the reaction to the list and increment the rID counter
                push!(reaction_list, forward_rxn)
                rID += 1
                                
                #### If the molecule is stable the reverse reaction will be VERY Slow (maybe 0?)
                if f1_grow in stable_molecules
                    backward_k = stable_backward
                else
                    backward_k = 1.0
                end
                #### Build the reverse reaction based on stability and add it to the list
                backward_rxn = Reaction(rID, [p1], [1], [m1,monomer_id], [1,1],[],[], "STD", backward_k)
                push!(reaction_list, backward_rxn)
                rID += 1 
                
                f1_cat = "Mo36*"*f1
                forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1_cat), count_mo_num("Mo1"))
                f1_grow_cat = "Mo36*"*f1_grow
                if (f1_cat in all_molecules) && (f1_grow_cat in all_molecules)
                    m1_cat = molecule_dict[f1_cat]
                    p1_cat = molecule_dict[f1_grow_cat]
                    cat_forward_rxn = Reaction(rID, [monomer_id,m1_cat], [1,1], [p1_cat], [1],[],[], "STD", forward_k)
                    push!(reaction_list, cat_forward_rxn)
                    rID += 1
                    cat_backward_rxn = Reaction(rID, [p1_cat],[1],[monomer_id,m1_cat], [1,1],[],[], "STD", 1.0/mo36_enhance)
                    push!(reaction_list, cat_backward_rxn)
                    rID += 1
                end
                
            end
            
        end
        
        for f2 in penta_fragments
            # Get the molecules id numbers
            if f2 != "Mo1"
                m2 = molecule_dict[f2]

                ## Set stability of (Mo6+2) vs Mo6
                # If the two molecules are not a set of reactants that have already been handled (This prevents you from making two reactions for A+B -> C and B+A->C)
                if !(Set([m1,m2]) in reactants_set)
                    ## Check if they're pentamer fragements that can be combined
                    if (!occursin(")_(", f1) && !occursin(")_(", f2) && all([f1,f2] .!= "Mo36" ))
                        ### Get their product
                        extra = 0
                        prod_f, extra = merge_penta_fragments(f1,f2)
                        ### If that product is allowed to be formed it is in all molecules
                        if (prod_f in all_molecules && !check_if_stackable(f1,f2))
                            #### Figure out what the product is and what it's formation reaction coefficient is
                            p1 = molecule_dict[prod_f]
                            #println("Reaction Number: ", rID, "     Reactants: ", [f1, f2], "   Product: ", prod_f)
                            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num(f2))
                            if prod_f in stable_molecules
                                backward_k = stable_backward
                            else
                                backward_k = 1.0
                            end
                            
                            if f1 != f2
                                if extra == 0
                                    forward_rxn = Reaction(rID, [m1,m2], [1,1], [p1], [1],[],[], "STD", forward_k)
                                    #### Addthe reaction to the list and increment the rID counter
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                  
                                    #### Build the reverse reaction based on stability and add it to the list
                                    backward_rxn = Reaction(rID, [p1], [1], [m1,m2], [1,1],[],[], "STD", backward_k)
                                    push!(reaction_list, backward_rxn)
                                    rID += 1
                                else
                                    forward_rxn = Reaction(rID, [m1,m2], [1,1], [p1, monomer_id], [1, extra],[],[], "STD", forward_k)
                                    #### Addthe reaction to the list and increment the rID counter
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    # #### Build the reverse reaction based on stability and add it to the list
                                    # backward_rxn = Reaction(rID, [p1, monomer_id], [1, extra], [m1,m2], [1,1],[],[], "STD", backward_k)
                                    # push!(reaction_list, backward_rxn)
                                    # rID += 1
                                end
                                
                            elseif f1 == f2
                                if extra == 0
                                    forward_rxn = Reaction(rID, [m1], [2], [p1], [1],[],[], "STD", forward_k)
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    #### Build the reverse reaction based on stability and add it to the list
                                    backward_rxn = Reaction(rID, [p1], [1], [m1], [2],[],[], "STD", backward_k)
                                    push!(reaction_list, backward_rxn)
                                    rID += 1
                                else
                                    forward_rxn = Reaction(rID, [m1], [2], [p1, monomer_id], [1, extra],[],[], "STD", forward_k)
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    #### Build the reverse reaction based on stability and add it to the list
                                    # backward_rxn = Reaction(rID, [p1, monomer_id], [1, extra], [m1], [2],[],[], "STD", backward_k)
                                    # push!(reaction_list, backward_rxn)
                                    # rID += 1
                                end
                            end
                            
                        elseif check_if_stackable(f1,f2)
                                extra = 0
                                prod_f,extra = stack(f1,f2)
                                p1 = molecule_dict[prod_f]
                                if prod_f in stable_molecules
                                        backward_k = stable_backward
                                else
                                    backward_k = 1.0
                                end
                                #println("Reaction Number: ", rID, "     Reactants: ", [f1, f2], "   Product: ", prod_f)
                                forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num(f2))
                                if f1 != f2
                                    if extra == 0
                                        forward_rxn = Reaction(rID, [m1,m2], [1,1], [p1], [1],[],[], "STD", forward_k)
                                        #### Add the reaction to the list and increment the rID counter
                                        push!(reaction_list, forward_rxn)
                                        rID += 1 
                                        #### Build the reverse reaction based on stability and add it to the list
                                        backward_rxn = Reaction(rID, [p1], [1], [m1,m2], [1,1],[],[], "STD", backward_k)
                                        push!(reaction_list, backward_rxn)
                                        rID += 1
                                    else
                                        forward_rxn = Reaction(rID, [m1,m2], [1,1], [p1, monomer_id], [1, extra],[],[], "STD", forward_k)
                                        push!(reaction_list, forward_rxn)
                                        rID += 1 
                                        #### Build the reverse reaction based on stability and add it to the list
                                        # backward_rxn = Reaction(rID, [p1, monomer_id], [1, extra], [m1,m2], [1,1],[],[], "STD", backward_k)
                                        # push!(reaction_list, backward_rxn)
                                        # rID += 1
                                    end
                                    
                                elseif f1 == f2
                                    if extra == 0
                                        forward_rxn = Reaction(rID, [m1], [2], [p1], [1],[],[], "STD", forward_k)
                                        #### Add the reaction to the list and increment the rID counter
                                        push!(reaction_list, forward_rxn)
                                        rID += 1 
                                        #### Build the reverse reaction based on stability and add it to the list
                                        backward_rxn = Reaction(rID, [p1], [1], [m1], [2],[],[], "STD", backward_k)
                                        push!(reaction_list, backward_rxn)
                                        rID += 1
                                    else
                                        forward_rxn = Reaction(rID, [m1], [2], [p1, monomer_id], [1, extra],[],[], "STD", forward_k)
                                        push!(reaction_list, forward_rxn)
                                        rID += 1 
                                        # #### Build the reverse reaction based on stability and add it to the list
                                        # backward_rxn = Reaction(rID, [p1, monomer_id], [1, extra], [m1], [2],[],[], "STD", backward_k)
                                        # push!(reaction_list, backward_rxn)
                                        # rID += 1
                                    end
                                end
                        end
                        #### Record which reactants you merged so you don't double count
                        push!(reactants_set, Set([m1,m2]))
                    ## Check if they're two stacks that can be merged to form Mo36
                    elseif ( occursin(")_(",f1) && occursin(")_(", f2) )
                        ### If the two stacks can be merged and their not already recorded
                       
                        merge_able, extra = can_merge_stacks(f1,f2)
                        if merge_able && !(Set([m1,m2]) in reactants_set)
                            #### Merged stacks always make Mo36
                            prod_f = "Mo36"
                            backward_k = stable_backward
                            #println("Reaction Number: ", rID, "     Reactants: ", [f1, f2], "   Product: ", prod_f)
                            #### Make reaction to form Mo36, stable, back rxn is almost non-existent, forward reaction is fast
                            p1 = molecule_dict[prod_f]
                            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num(f2))
                            if f1 != f2
                                if extra == 0
                                    forward_rxn = Reaction(rID, [m1,m2], [1,1], [p1], [1],[],[], "STD", forward_k)
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    backward_rxn = Reaction(rID, [p1], [1], [m1,m2], [1,1],[],[], "STD", backward_k)
                                    push!(reaction_list, backward_rxn)
                                    rID += 1
                                else
                                    forward_rxn = Reaction(rID, [m1,m2], [1,1], [p1, monomer_id], [1, extra],[],[], "STD", forward_k)
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    # backward_rxn = Reaction(rID, [p1, monomer_id], [1, extra], [m1,m2], [1,1],[],[], "STD", backward_k)
                                    # push!(reaction_list, backward_rxn)
                                    # rID += 1
                                end
                                    
                            elseif f1 ==f2
                                if extra == 0 
                                    forward_rxn = Reaction(rID, [m1], [2], [p1], [1],[],[], "STD", forward_k)
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    backward_rxn = Reaction(rID, [p1], [1], [m1], [2],[],[], "STD", backward_k)
                                    push!(reaction_list, backward_rxn)
                                    rID += 1
                                else
                                    forward_rxn = Reaction(rID, [m1], [2], [p1, monomer_id], [1, extra],[],[], "STD", forward_k)
                                    push!(reaction_list, forward_rxn)
                                    rID += 1
                                    # backward_rxn = Reaction(rID, [p1, monomer_id], [1, extra], [m1], [2],[],[], "STD", backward_k)
                                    # push!(reaction_list, backward_rxn)
                                    # rID +=1
                                end
                            end
                            push!(reactants_set, Set([m1,m2]))
                        end

                    ## Check if one is a stack that can add a monomer
                    elseif f1== "Mo1" && occursin(")_(",f2)
                        if can_add_monomer_to_stack(f2)
                            possible_products = add_monomer_to_stack_all(f2)
                            for p in possible_products
                                if p in stable_molecules
                                    backward_k = stable_backward
                                else
                                    backward_k = 1.0
                                end
                                ## Figure out what the product is an how stable it is
                                p1 = molecule_dict[p]
                                #println("Reaction Number: ", rID, "     Reactants: ", [f1, f2], "   Product: ", p)
                                forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f2), count_mo_num("Mo1"))
                                forward_rxn = Reaction(rID, [m2,monomer_id], [1,1], [p1], [1],[],[], "STD", forward_k)
                                push!(reaction_list, forward_rxn)
                                rID += 1  
                                backward_rxn = Reaction(rID, [p1], [1], [m2,monomer_id], [1,1],[],[], "STD", backward_k)
                                push!(reaction_list, backward_rxn)
                                rID += 1
                            end
                            push!(reactants_set, Set([m1,m2]))
                        end
                    end
                end
            end
        end
    end
    
    # Next make all the reactions to assemble the wheel
    for f1 in wheel_fragments
        m1 = molecule_dict[f1]
        # Add a pentamer to this fragment
        prod_f = add_pent_wheel(f1)
        # Check if that product is possible and if the reaction hasn't be recorded yet
        if prod_f in all_molecules && !(Set([m1,penta_id]) in reactants_set)
            ## Whats the product and record the reaction
            p1 = molecule_dict[prod_f]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo6"], "   Product: ", prod_f)
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num("Mo6"))
            forward_rxn = Reaction(rID, [m1,penta_id], [1,1], [p1], [1],[],[], "STD", calculate_wheel_stability(forward_k, wheel_enhancement_multiplier, prod_f))
            push!(reaction_list, forward_rxn)
            rID += 1

            if prod_f in stable_molecules
                backward_k = stable_backward
            else
                backward_k = 1.0
            end
            
            backward_rxn = Reaction(rID, [p1], [1], [m1,penta_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
            
            push!(reactants_set, Set([m1, penta_id]))
        end
        # A a monomer to the wheel and see if its a possible product that hasn't been recorded yet
        prod_f = add_monomer_wheel(f1)
        if prod_f in all_molecules && !(Set([m1,monomer_id]) in reactants_set)
            
            p1 = molecule_dict[prod_f]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo1"], "   Product: ", prod_f)
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num("Mo1"))
            forward_rxn = Reaction(rID, [m1,monomer_id], [1,1], [p1], [1],[],[], "STD", calculate_wheel_stability(forward_k, wheel_enhancement_multiplier, prod_f))
            push!(reaction_list, forward_rxn)
            rID += 1

            if prod_f in stable_molecules
                backward_k = stable_backward
            else
                backward_k = 1.0
            end
            backward_rxn = Reaction(rID, [p1], [1], [m1,monomer_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
            
            push!(reactants_set, Set([m1, monomer_id]))
        end
        # Add a dimer to the wheel fragment and see if it's possible and whether it's been recorded 
        prod_f = add_dimer_wheel(f1)
        if prod_f in all_molecules && !(Set([m1,dimer_id]) in reactants_set)
            p1 = molecule_dict[prod_f]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo2"], "   Product: ", prod_f)
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num("Mo2"))
            forward_rxn = Reaction(rID, [m1,dimer_id], [1,1], [p1], [1],[],[], "STD", calculate_wheel_stability(forward_k, wheel_enhancement_multiplier, prod_f))
            push!(reaction_list, forward_rxn)
            rID += 1

            if prod_f in stable_molecules
                backward_k = stable_backward
            else
                backward_k = 1.0
            end
            
            backward_rxn = Reaction(rID, [p1], [1], [m1,dimer_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
            
            push!(reactants_set, Set([m1, dimer_id]))
        end
                
    end
    
    # Repeat for the ball fragments, remember balls are built using edgeMo2 dimers
    for f1 in ball_fragments
        m1 = molecule_dict[f1]
        prod_f = add_pentamer_ball(f1)
        if prod_f in all_molecules && !(Set([m1,penta_id]) in reactants_set)
            p1 = molecule_dict[prod_f]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "Mo6"], "   Product: ", prod_f)
            k_ball_f = calculate_ball_addition_rate(k_f, ball_growth_multipler, f1, "Mo6")
            forward_k = bimolecular_coef(k_ball_f,T,R,volume,count_mo_num(f1), count_mo_num("Mo6"))
            forward_rxn = Reaction(rID, [m1,penta_id], [1,1], [p1], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1

            if prod_f in stable_molecules
                backward_k = stable_backward
            else
                backward_k = 1.0
            end
            backward_rxn = Reaction(rID, [p1], [1], [m1,penta_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
            
            push!(reactants_set, Set([m1, penta_id]))
        end

        prod_f = add_dimer_ball(f1)
        if prod_f in all_molecules && !(Set([m1,edge_dimer_id]) in reactants_set)
            p1 = molecule_dict[prod_f]
            #println("Reaction Number: ", rID, "     Reactants: ", [f1, "edgeMo2"], "   Product: ", prod_f)
            k_ball_f = calculate_ball_addition_rate(k_f, ball_growth_multipler, f1, "edgeMo2")
            forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num(f1), count_mo_num("edgeMo2"))
            forward_rxn = Reaction(rID, [m1,edge_dimer_id], [1,1], [p1], [1],[],[], "STD", forward_k)
            push!(reaction_list, forward_rxn)
            rID += 1

            if prod_f in stable_molecules
                backward_k = stable_backward
            else
                backward_k = 1.0
            end
            backward_rxn = Reaction(rID, [p1], [1], [m1,edge_dimer_id], [1,1],[],[], "STD", backward_k)
            push!(reaction_list, backward_rxn)
            rID += 1
            push!(reactants_set, Set([m1, edge_dimer_id]))
        end           
    end
    
    #### Build the association of Mo36 and a pentamer
    prod_f = "Mo36*(Mo6)1_(Mo2)0_(Mo1)0"
    p = molecule_dict[prod_f]
    forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num("Mo6"), count_mo_num("Mo36"))
    forward_reaction = Reaction(rID, [mo36_id, penta_id], [1,1], [p], [1], [], [], "STD", forward_k)
    push!(reaction_list, forward_reaction)
    rID += 1
    backward_reaction = Reaction(rID, [p], [1], [mo36_id, penta_id], [1,1], [], [], "STD", 1.0)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    ### Build the association of pentamer and edgeMo2
    prod_f = build_ball_frag(1,1)
    p = molecule_dict[prod_f]
    forward_k = bimolecular_coef(k_f,T,R,volume,count_mo_num("Mo6"), count_mo_num("edgeMo2"))
    forward_reaction = Reaction(rID, [edge_dimer_id, penta_id], [1,1], [p], [1], [], [], "STD", forward_k)
    push!(reaction_list, forward_reaction)
    rID += 1
    backward_reaction = Reaction(rID, [p], [1], [edge_dimer_id, penta_id], [1,1], [], [], "STD", 1.0)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    #### Build the dissociation of Mo36 from the wheel this reaction goes one way
    associated_id =  molecule_dict["Mo36*(Mo6)14_(Mo2)28_(Mo1)14"]
    
    wheel_id = molecule_dict["(Mo6)14_(Mo2)28_(Mo1)14"]
    backward_rxn = Reaction(rID, [associated_id], [1], [mo36_id, wheel_id], [1,1],[],[], "STD", 100.0)
    push!(reaction_list, backward_rxn)
    rID += 1
    associated_id =  molecule_dict["Mo36*Mo6"]
    
    
    ### Build the association Mo6 
    #### Build the dissociation of Mo6 from the Mo36 this reaction goes one way
    penta_id = molecule_dict["Mo6"]
    backward_rxn = Reaction(rID, [associated_id], [1], [mo36_id, penta_id], [1,1],[],[], "STD", 100.0)
    push!(reaction_list, backward_rxn)
    rID += 1
    
    parameters = Dict(
    :k_f => k_f,
    :stable_backward=> stable_backward,
    :dimerization_ratio => dimerization_ratio,
    :mo36_enhance => mo36_enhance,
    :wheel_enhancement_multiplier => wheel_enhancement_multiplier,
    :ball_enhance => ball_growth_multipler,
    :volume => volume,
    :T => T,
    :R => R
    )
    
    MoBlue_CRS = CRS(all_molecules, molecule_dict, reaction_list, parameters)
    
    
    return MoBlue_CRS
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
    elseif (core <5 && additions >=3)
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
    extra = 0 # extra monomers
    m_new = "NULL"
    
    mo1,add1,hanger1 = count_components_penta(m1)
    mo2,add2,hanger2 = count_components_penta(m2)
    mo_sum = sum([mo1, add1, hanger1, mo2, add2, hanger2])
    if mo_sum >= 10
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