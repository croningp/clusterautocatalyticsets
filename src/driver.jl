@everywhere using Random
@everywhere using DataFrames
@everywhere using JSON
@everywhere using CSV
@everywhere using Combinatorics
@everywhere include("Gillespie.jl")

using Dates
function runif(min::Float64,max::Float64, n::Int64 = 1)
    r = (max - min).*rand(n) .+ min
    return r
end

function logunif(min::Float64,max::Float64, n::Int64 = 1)
    exps = (max - min).*rand(n) .+ min
    r = 10.0 .^ exps
    return r
end

function Sample_MoBlueParameters(total_time::Float64, n_sims::Int64, seed::Int64)
    Random.seed!(seed)
    seeds = rand(seed:(seed + n_sims), n_sims)
    total_mass = 1000000
    # You are trying a very extreme set of parameters, very low "stable_backward" and higher k_f, to see if this traps the mo36 structure 
    sleep_time = logunif(-1.0, 1.0, n_sims)
    k_f = 0.3*1e-3 #logunif(-6.0,0.0, n_sims)
    k_d_stable = 0.01# logunif(-6.0,0.0, n_sims)
    dimerization_ratio = logunif(-3.0,2.0, n_sims)
    mo6_enhance = 0.0#logunif(-1.0, 2.0, n_sims)
    mo36_enhance = 100.0#logunif(0.0, 3.0, n_sims)
    #nano_stuct_enhance =0.0#
    nano_stuct_enhance = 100.0#logunif(1.0, 3.0, n_sims)
    @sync @distributed for i=1:n_sims
        sleep(sleep_time[i])
        #i_seed = Int64(floor(10.0^(10^sleep_time[i])))
        Run_MoBlue(seeds[i], total_mass, total_time, k_f, k_d_stable, dimerization_ratio[i], mo6_enhance, mo36_enhance, nano_stuct_enhance, nano_stuct_enhance)
    end
    
    
end