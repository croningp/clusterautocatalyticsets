@everywhere using Random
@everywhere using DataFrames
@everywhere using JSON
@everywhere using CSV
@everywhere using Combinatorics
@everywhere include("Gillespie.jl")


function logunif(min::Float64,max::Float64, n::Int64 = 1)
    exps = (max - min).*rand(n) .+ min
    r = 10.0 .^ exps
    return r
end

function Sample_MoBlueParameters(total_mass::Int64, total_time::Float64, n_sims::Int64)
    
    
    # You are trying a very extreme set of parameters, very low "stable_backward" and higher k_f, to see if this traps the mo36 structure 
    sleep_time = logunif(-2.0, 0.0, n_sims)
    k_f =1E-4 #logunif(-6.0, -1.0)
    stable_backward = 1E-8#logunif(-6.0,0.0, n_sims)
    dimerization_ratio = logunif(-0.5, 1.0, n_sims)
    mo36_enhance = 0.1#rand([1.0, 10.0, 100.0], n_sims)#
    ball_growth_multiplier =0.0 #logunif(0.0, 6.0)
    wheel_growth_multiplier =0.0
    
    @sync @distributed for i=1:n_sims
        sleep(sleep_time[i])
        Run_MoBlue(total_mass, total_time, k_f, stable_backward, dimerization_ratio[i], mo36_enhance, wheel_growth_multiplier, ball_growth_multiplier)
    end
    
    
end