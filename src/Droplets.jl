#############################################
############ Droplet Functions ##############
#############################################

mutable struct Droplet # Individual droplets
    ID::Int64
    xpos::Int64
    ypos::Int64
    concentrations::Array{Int64,1}
    volume::Float64
    CRS::CRS
end


function initialize_droplets(num_drops::Int64, lattice_L::Int64, max_L::Int64, kf::Float64, kb::Float64, init_monomers::Int64)
    
    volume = 1.0/float(num_drops)
    lattice = zeros(Int64, lattice_L,lattice_L)
    x,y = 1:lattice_L, 1:lattice_L
    xypairs = [repeat(x, inner=[size(y,1)]) repeat(y, outer=[size(x,1)])]
    my_pairs = rand(1:size(xypairs)[1], num_drops)
    xs = xypairs[my_pairs,1]
    ys = xypairs[my_pairs,2]
    drops = Array{Droplet}(undef, num_drops)
    
    original_CRS = generate_binary_ligation_CRS(max_L, kf, kb, volume, 300.0, 8.3144598, 100, 100)
    num_molecules = length(original_CRS.molecule_list)
    Aindex = original_CRS.molecule_dict["A"]
    Bindex = original_CRS.molecule_dict["B"]
    
    
    for i in 1:num_drops
        init_A = floor(Random.rand()*init_monomers)
        init_B = init_monomers - init_A
        init_conc = zeros(num_molecules)
        init_conc[Aindex] = init_A
        init_conc[Bindex] = init_B
        lattice[xs[i], ys[i]] = i
        drops[i] = Droplet(i, xs[i], ys[i], init_conc, volume, deepcopy(original_CRS))
        
    end
    
    return drops,lattice
    
end

function time_evolve_droplets(lattice::Array{Int64,2}, drops::Array{Droplet,1}, time, delta_t, split_p)
    n = Int(time/delta_t)
    L = size(lattice)[1]
    t= 0.0
    
    conc_DF = DataFrame(droplet =Int64[], t = Float64[],conc = Array{Int64,1}[])
    drop_edge_DF = DataFrame(t= Float64[], source = Int64[], target = Int64[])
    ### Time evolution loop
    for i in 1:n ### For all the droplets update their internal chemistry
        println(i)
        #Threads.@threads
        for m in 1:length(drops)
            drops[m]= evolve_droplet_chemistry(drops[m], t, t+delta_t)
        end
        
        
        for d1 in drops ### For all the droplets move or split
            push!(conc_DF, [d1.ID, i*delta_t, d1.concentrations])
            new_x,new_y = pick_next(d1.xpos, d1.ypos, L) ## Get the next position
            
            if (lattice[new_x,new_y] != 0.0) ## See if it's occupied 
                #### There's definitely a bug here. 
                
                d2 = drops[1]
                
                # If it's occupied find the occupant and combine
                found = false
                j = 1
                while found == false
                   if drops[j].xpos == new_x && drops[j].ypos == new_y && drops[j].ID != d1.ID
                        d2= drops[j]
                        found = true
                    else
                       j += 1
                    end
                end
                
                # Combine, remove old drops, and add new one to drops array
                new_drop,lattice = combine_droplets(d1,d2,lattice)
                push!(drop_edge_DF, [i*delta_t, d1.ID, new_drop.ID],[i*delta_t, d2.ID, new_drop.ID] )
                let d1=d1, d2=d2 # Hassam magic to solve the a type instability Bug 15276
                    filter!(e-> e != d1 && e !=d2, drops)
                end
                pushfirst!(drops, new_drop)     
                
            else
                ### If it's not occupied 
                ## Check for split event
                if Random.rand() <= delta_t*split_p ## Split
                    
                    new1,new2,lattice = split_droplet(d1, lattice, new_x, new_y)
                    push!(drop_edge_DF, [i*delta_t, d1.ID, new1.ID],[i*delta_t, d1.ID, new2.ID] )
                    let d1=d1 # Hassam magic to solve the a type instability Bug 15276
                        filter!(e-> e ∉ [d1], drops)
                    end
                    pushfirst!(drops, new1, new2)
                else
                     ## Move
                    
                    lattice[new_x,new_y] = d1.ID
                    lattice[d1.xpos, d1.ypos] = 0.0
                    d1.xpos = new_x
                    d1.ypos = new_y
                    push!(drop_edge_DF, [i*delta_t, d1.ID, d1.ID])
                end 
            end ### If it's not occupied 
        end ### Drop loop
        
    end ### Time Loop

    return conc_DF, drop_edge_DF
end



function pick_next(x::Int64,y::Int64, L::Int64)
    r = Random.rand(1:4)
    if r == 1 # Up
        y += 1
    elseif r ==2 # Down
        y -= 1
    elseif r == 3 # Left
        x -=1
    elseif r == 4 # Right
        x +=1
    
    end
    x = mod1(x, L)
    y = mod1(y, L)
    
   return x,y 
end

function split_concentrations(c)
    c1 = zero(c)
    c2 = zero(c)
    for n in 1:length(c)
        a= Random.rand(0:c[n])
        c1[n] = a
        c2[n] = c[n] -a
    end
    return c1,c2
end



function split_droplet(drop::Droplet, lattice::Array{Int64,2}, new_x::Int64, new_y::Int64)
    maxID, index = findmax(lattice)
    
    c1, c2 = split_concentrations(drop.concentrations)
    CRS1 = copy_CRS(drop.CRS, 0.5)
    CRS2 = copy_CRS(drop.CRS, 0.5)
    
    lattice[drop.xpos, drop.ypos] = maxID + 1
    drop1 = Droplet( (maxID+1), drop.xpos, drop.ypos, c1, 0.5*drop.volume, CRS1 )
    
    lattice[new_x,new_y] = maxID +2 
    drop2 = Droplet( (maxID+2), new_x, new_y, c2, 0.5*drop.volume, CRS2)
    return drop1, drop2, lattice
end

function combine_droplets(drop1,drop2, lattice)
    maxID, idx = findmax(lattice)
    new_volume_ratio = (drop1.volume + drop2.volume)/ drop1.volume
    new_CRS = copy_CRS(drop1.CRS, new_volume_ratio)
    new_c = drop1.concentrations + drop2.concentrations
    
    lattice[drop2.xpos, drop2.ypos] = 0
    drop = Droplet(maxID+1, drop1.xpos, drop1.ypos, new_c, new_volume_ratio*drop1.volume, new_CRS)
    lattice[drop1.xpos, drop1.ypos] = maxID+1
    return drop, lattice 
end



function evolve_droplet_chemistry(drop, t, next_t)
    #Initialize Propensities
    concentrations = drop.concentrations[:]
    propensities = CalculatePropensities(concentrations, drop.CRS)
    ####### Main LOOP #######
    while t < next_t
        # Pick Reaction
        rID = PickReactionID(propensities)
        # Execute Reaction
        concentrations = ExecuteReaction(concentrations, drop.CRS, rID)
        # Calculate Propensities
        propensities = CalculatePropensities(concentrations, drop.CRS)
        
        #Update Time
        Ap_tot = sum(propensities)
        t -= (log(rand())/Ap_tot) 
    end
    drop.concentrations = concentrations[:]
   return drop 
end
