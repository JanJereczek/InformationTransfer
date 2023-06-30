using DrWatson
@quickactivate "InformationTransfer"

using LinearAlgebra, Random, Statistics
using NCDatasets, JLD2, DelimitedFiles
using CairoMakie, GeoMakie
using TSVD, Interpolations

# Here you may include files from the source directory
include(srcdir("core.jl"))
include(srcdir("utils.jl"))
include(srcdir("plotting.jl"))
