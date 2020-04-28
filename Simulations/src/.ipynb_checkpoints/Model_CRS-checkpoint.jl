#using Iterators

#############################################
###### Structs ##############################
#############################################
struct Reaction ## Reaction Object
	index::Int64 # List index of reaction
	reactants::Array{Int64,1} # List of reactant indices
	react_coef::Array{Int64,1} # List of reactant coefficients
	products::Array{Int64,1} # List of product indices
	prod_coef::Array{Int64,1} # List of product coefficients
	catalysts::Array{Int64,1} # Catalyst indicies
	cat_coef::Array{Float64,1} # Catalysts coefficients
	propensity::String # Propensity identifier
	rate_constant::Float64 # Rate constant
end

struct CRS # Chemical reaction system object
	molecule_list::Array{String,1} # List of molcule strings
	molecule_dict::Dict{String, Int64} # Dict maps molecule strings to index
	reaction_list::Array{Reaction,1} # List of reaction objects
	parameters::Dict{Any,Any} # Every other parameter
end


#############################################
###### CRS functions ########################
#############################################


function generate_all_binary_seqs(max_L)
    ### Generate all possible sequences
    all_seqs = []
    for j in 1:max_L
        comps = Array{String}(undef,j + 1)
        for l in 0:j
            comps[l+1] = "B"^l * "A"^(j-l)
        end
        
        for c in comps
            p = unique(collect(permutations(c)))
            all_seqs = vcat(p,all_seqs)
        end         
    end
    seqs = [join(s) for s in all_seqs]    
    sort!(seqs, by = length, rev= false)
    return seqs
end

### Count occurance of t in s Function ##
function countall(s, t)
    n = length(collect(eachmatch(Regex(t), s)))
    return n
end

function all_splits(seq)
    l = length(seq)
    left = Array{String}(undef,(l-1))
    right =Array{String}(undef,(l-1))
    
    for i in 1:(l-1)
        left[i] = seq[1:i]
        right[i] = seq[(i+1):end] 
    end
    return left,right
end

function calculate_reduced_mass_binary(m1_string, m2_string, mA,mB)
    m1_mass = countall(m1_string,"A")*mA + countall(m1_string,"B")*mB
    m2_mass = countall(m2_string,"A")*mA + countall(m2_string,"B")*mB
    reduce = reduced_mass(m1_mass,m2_mass)
    return reduce
end

function reduced_mass(m1,m2)
    return (m1*m2)/(m1+m2)
end

function bimolecular_coef(k, T, R, V,m1, m2)
    m = reduced_mass(m1,m2)
    return (k/V)*sqrt(R*T/m)
end 

function generate_binary_ligation_CRS(max_L, kf, kb, volume::Float64 = 1.0, T::Float64 = 300, R::Float64= 8.3144598, mA::Int64= 100, mB::Int64 = 100)
    seqs = generate_all_binary_seqs(max_L)
    ## Generate Molecule List and Dictionary
    molecule_dict = Dict{String,Int64}()
    for (n, m) in enumerate(seqs)
        molecule_dict[m] = n
    end
    
    ## Generate Reaction List
    Aindex = molecule_dict["A"]
    Bindex = molecule_dict["B"]
    
    ### Generate all forward and backward reactions
    n_rxn = sum([(i-1)*2^(i +1) for i in 1:max_L  ])
    rxn_list = Array{Reaction}(undef,n_rxn )
    rID = 1
    for s in seqs
        ## Get all splits of the sequence 
        left_side, right_side = all_splits(s)
        for i in 1:length(left_side)
            ### Forward Reaction
            ## Check if both sides of sequence are identical 
            if left_side[i] != right_side[i]
                reactants = [ molecule_dict[ left_side[i] ], molecule_dict[ right_side[i] ] ]
                reactant_coef = [1,1]
            else
                reactants = [molecule_dict[ left_side[i] ] ]
                reactant_coef = [2]
            end
            
            products = [ molecule_dict[s] ]
            product_coef = [1]
            
            m_red = calculate_reduced_mass_binary(left_side[i], right_side[i],mA, mB)
            k = bimolecular_coef(kf, T, R, volume, m_red)
            
            forward_rxn = Reaction(rID, reactants, reactant_coef, products, product_coef, Array{Int64}(undef,0), Array{Float64}(undef,0), "STD", k)
            rxn_list[rID] = forward_rxn
            rID += 1
            
            ## Backward Reaction
            backward_rxn = Reaction(rID, products, product_coef,reactants, reactant_coef, Array{Int64}(undef,0), Array{Float64}(undef,0), "STD", kb)
            rxn_list[rID] = backward_rxn
            rID += 1
            
        end
            
    end  
    
    parameters = Dict(:max_L => max_L, :kf => kf, :kb => kb, :volume => volume, :T => T, :R => R, :mA=> mA, :mB => mB)
    return(CRS(seqs, molecule_dict, rxn_list, parameters))
end

function copy_CRS(original_CRS, volume_multiplier)
    reaction_list = similar(original_CRS.reaction_list)
    for rID in 1:length(reaction_list)
        rxn = original_CRS.reaction_list[rID]
        if sum(rxn.react_coef) ==1
            reaction_list[rID] = deepcopy(original_CRS.reaction_list[rID])
        elseif sum(rxn.react_coef) == 2
            new_rxn = Reaction(rxn.index, rxn.reactants, rxn.react_coef, rxn.products, rxn.prod_coef, rxn.catalysts, rxn.cat_coef, rxn.propensity, (1.0/volume_multiplier)*rxn.rate_constant)
            reaction_list[rID] = new_rxn
        end
       
    end
    new_parameters = original_CRS.parameters
    new_parameters[:volume] *= volume_multiplier
    
    new_CRS = CRS(original_CRS.molecule_list, original_CRS.molecule_dict, reaction_list, new_parameters )
    return new_CRS 
end

function generate_tuples(max, sum)
    pairs = collect(Iterators.product(1:max, 1:max))
    my_pairs = []
    for pa in pairs
        if sum(pa) <= sum
            if !(pa in my_pairs) && !(reverse(pa) in my_pairs)
                push!(my_pairs,pa)
            end
        end
    end
    return my_pairs
end

