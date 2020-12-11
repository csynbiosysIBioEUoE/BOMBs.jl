module BOMBS

# Load dependencies
# using BayesianOptimization
# using BlackBoxOptim
using CSV
# using Calculus
# using CmdStan
# using DataFrames
using Dates
using DiffEqBase
using DifferentialEquations
# using Distributed
# using Distributions
# using GaussianProcesses
using JLD
# using LinearAlgebra
using ODEInterfaceDiffEq
using OrdinaryDiffEq
using Plots
# using Random
# using Roots
# # using SharedArrays
# using Statistics
# using StatsBase
# using StatsPlots
using Sundials

# laod logo
include("logoBOMBS.jl")
export printLogo
export versionBOMBS

# load info function
include("InfoStructs.jl")
export infoAll

# Load section to generate the model scripts
include("ModelGen.jl")
export defModStruct
export checkStruct
export GenerateModel

include("ModelSim.jl")
export simulateODEs
export defSimulStruct
export checkStructSimul
export fileStructInfo
export defSimulStructFiles
export extractSimulCSV
export plotSimsODE


end
