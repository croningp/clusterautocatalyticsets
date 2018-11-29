# ChemEvolve.jl
Dr Cole Mathis, Dr Sara Walker, ELife@ASU & CroninLab
September 2018

Julia Implementation of Gillespie Algorithm with flexible tools for dealing 
with stochastic reaction systems. 





## Introduction

This software package will enable the flexible implementation of the 
Gillespie algorithm in Julia, along with analysis tools for 
characterizing complex chemical reaction systems. 

Many different chemical models will be implemented, including:
Ligation and Polymerization Chemistry
Sequence Based Replication Chemistry
Reflexively Autocatalytic and Food Generated Sets
The Azoarcus Recombinase System
and other random autocatalytic networks

## Kinetics

Bimolecular reactions occur with a rate constant given by
```math
k_r = \frac{c_r}{V} \sqrt{ \frac{R T}{m_ {ij}} }
```

where $`c_r`$ corresponds to some structural and kinetic factors which are 
unknown, $`T`$ is the temperature, $`m_{ij}`$ is the reduced mass in daltons.
of the two reactants, $`V`$ is the reaction volume and $`R`$ is the ideal gas \
constant.

