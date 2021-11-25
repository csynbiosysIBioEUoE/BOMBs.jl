module BOMBs

# Load dependencies
using BayesianOptimization
using BlackBoxOptim
using CSV
using Calculus
using CmdStan
using DataFrames
using Dates
using DiffEqBase
using DifferentialEquations
using Distributed
using Distributions
using GaussianProcesses
using GaussianMixtures
using JLD
using LinearAlgebra
# using MCMCChains
using ODEInterfaceDiffEq
using OrdinaryDiffEq
using Plots
using Random
using Roots
using ScikitLearn.GridSearch
using ScikitLearnBase
using SharedArrays
# using Statistics
using StatsBase
using StatsPlots
using Sundials
# using Turing

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

include("ModelStanInfer.jl")
export defBayInfStruct
export defBayInfDataStruct
export defBayInfDataFromFilesStruct
export defBasicStanSettingsStruct
export convertBoundTo2
export convertBoundToReal
export fitPriorSamps
export fitPriorSampsMultiNorm
export checkStructBayInf
export checkStructBayInfData
export checkStructBayInfDataFiles
export checkStructBayInfStanSettings
export genStanInitDict
export reparamDictStan
export genStanModel
export restructureDataInference
export getStanInferenceElements
export saveStanResults
export runStanInference
export plotStanResults
export StanInfer

include("ModelTuringInfer.jl")
export defTurInfStruct
export fitPriorSampsMultiNormTuring
export checkStructTurInf
export genTuringModel
export getTuringInferenceElements
export saveTuringResults
export plotTuringResults
export TuringInfer

include("ModelEntropyTheta.jl")
export genSamplesPrior
export computeH
export computeHgain

include("ModelOEDSelection.jl")
export defODEModelSelectStruct
export checkStructOEDMS
export BhattacharyyaDist
export EuclideanDist
export genOptimMSFuncts
export plotOEDMSResults
export settingsBayesOpt
export mainOEDMS

include("ModelOEDCalibration.jl")
export defODEModelCalibrStruct
export checkStructOEDMC
export genOptimMCFuncts
export plotOEDMCResults
export settingsBayesOptMC
export mainOEDMC


end
