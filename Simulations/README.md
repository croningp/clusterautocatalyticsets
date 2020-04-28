# Gillespie Mo Blue Simulations
Dr Cole Mathis, Dr Harry Miras, Dr Lee Cronin

Feb 2019

Julia Implementation of Gillespie Algorithm simulation of Mo-Blue kinetic model


# Dependencies 
This code uses Julia (v1.0 or later) for simulations and R (v3.5 or later) for analysis
The required Julia packages are

- DataFrames
- CSV
- Combinatorics
- JSON
- Random

The required R packages are 

- ggplot2
- dplyr
- tidyr
(all of these, and much more, are installed with `install.packages('tidyverse')`)

The analysis is done using R via a Jupyter notebook which enables interactive development. More information on install JupyterLab can be found here
<https://jupyterlab.readthedocs.io/en/latest/>

More information on Julia, R, and tidyverse are availible at <https://julialang.org/>, <https://www.r-project.org/>, <https://www.tidyverse.org/>

The analysis is done using a Jupyter notebook with the IR kernal information for installation at <https://irkernel.github.io/installation/>


# Running Simulations
The Julia implementation allows the Simulation to be run with a single initial condition, or in parallelized batches. 

## Single Runs
In order to run a single instance read the main code file into Julia `include("Gillespie.jl")`

The function to call is named `Run_MoBlue(seed::Int64,Mo1_mass::Int64, total_time::Float64, k_f::Float64, k_d_stable::Float64, dimerization_ratio::Float64, mo6_enhance::Float64, mo36_enhance::Float64, mo154_enhance::Float64= 1.0,mo132_enhance::Float64 = 1.0)`

It has a large number of parameters, some model parameters and some simulation conditions. 

* seed: An integer which sets the random seed for the simulations. The stochastic nature of the simulation means that different random seeds generate different trajectories even for the same initial condition
* Mo1_mass: The total number of Molybdate molecules
* total_time: total simulation time in dimensionless units
* k_f: base forward reaction constant (see below for formula used to calculate rates)
* k_d_stable: the degradation rate constant for stable molecules (such as Mo2, Mo6, Mo36, Mo154, Mo132, etc...)
* dimerization_ratio: the ratio between edge-bonded Mo2 rate constant and corner bonded Mo2 formatio. The formation rate of edge-bonded Mo2 is proportional dimerization_ratio`*`k_f
* mo6_enhance: the degree to which Mo6 stablizes the formation of Mo6 molecules (autocatalytic small molecules)
* mo36_enhance: the degree to which Mo36 stablizes the formation of Mo6 molecules
* mo154_enhance: the degree to which mo154 fragments grow faster as they are completed
* mo132_enhance: the degree to which mo132 fragments grow faster as they are completed

Running the a single simulation will generate the timeseries data and save a csv into a folder called "/data/MoBlue_timeseries". The name of the file is a hashed-key. A seperate file maps the parameters run to the hashed filename. That file is called "data/MoBlue_Run_Parameters.csv". *This system may seem inefficient, however for paralellized batches of simulations this enables a much easier analysis work-flow*.

## Parallel Runs
Given the large number of parameters, it's often best to randomly sample the parameter space in order to characterize the dynamics of the system. To do that a Julia session needs to be started with multiple processors. If Julia isn't running launch it with the multiprocessor flag, for example `julia -p 10` will start a session with 10 processors. It's wise to pick the number of processors to be less than the maximum number of processors availible on your machine. For example, this software was developed on a PC with 40 real processing cores, so the number of processors used in Julia was typically set to 30. 

The file `driver.jl` contains functions to randomly sample parameters from various ranges and generate time series from those parameters. This allows us to quickly search the parameter space. To do this modify the `Sample_MoBlue(total_time, n_sims, seed)` function in the `driver.jl` script to vary the appropriate parameters. Some parameters naturally vary on logarithmic scales, the `logunif(min,max,n)` function will generate `n` values which are sampled uniformly on a log scale on the range `(10^min, 10^max)`. Once the function has been modified to sample the parameter ranges of interest the parallel batch can be submitted by loading the `driver.jl` file into Julia and calling the `Sample_MoBlue(total_time, n_sims, seed)` function. 

* total_time: the total simulation time of *each* simulation
* n_sims: the number of simulations (if this is greater than the number processors availible they will be submitted contineously)
* seed: random seed to control the initial conditions

