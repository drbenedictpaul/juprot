ENV["JULIA_PYTHONCALL_EXE"] = joinpath(ENV["HOME"], "plip_venv", "bin", "python3")
ENV["JULIA_CONDAEXE"] = ""
using Pkg
Pkg.activate(".")
include("src/JuProtGUI.jl")
using .JuProtGUI
println("JuProtGUI loaded successfully")
