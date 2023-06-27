using DrWatson
@quickactivate "InformationTransfer"

using LinearAlgebra, Random, Statistics
using NCDatasets, JLD2, DelimitedFiles
using CairoMakie, GeoMakie

# Here you may include files from the source directory
include(srcdir("core.jl"))
include(srcdir("utils.jl"))