module BOMBS

# Load dependencies
# using BayesianOptimization
using BlackBoxOptim
using CSV
# using Calculus
# using CmdStan
using DataFrames
using Dates
using DiffEqBase
using DifferentialEquations
using Distributed
using Distributions
# using GaussianProcesses
using JLD
using LinearAlgebra
using ODEInterfaceDiffEq
using OrdinaryDiffEq
using Plots
using Random
# using Roots
using SharedArrays
# using Statistics
using StatsBase
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

# Load section to simulate the model
include("ModelSim.jl")
export simulateODEs
export defSimulStruct
export checkStructSimul
export fileStructInfo
export defSimulStructFiles
export extractSimulCSV
export plotSimsODE

# Load the section to generate pseudo-data
include("ModelPDat.jl")
export GenPseudoDat
export plotPseudoDatODE
export defPseudoDatStruct
export checkStructPseudoDat
# export fileStructInfo
export defPseudoDatStructFiles
export extractPseudoDatCSV
export PDatCSVGen

include("ModelMLE.jl")
export defMLEStruct
export SimToMle
export checkStructMLE
export selectObsSim_te
export restructInputs_te
export UVloglike
export MVloglike
export plotMLEResults
export defCrossValMLEStruct
export checkStructCrossValMLE
export plotCrossValMLEResults
export CrossValMLE
export finishMLEres
export MLEtheta

end
