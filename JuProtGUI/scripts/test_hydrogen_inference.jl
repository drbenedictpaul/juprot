using BioStructures
using LinearAlgebra
using Statistics

# Include the hydrogen inference function
include("../lib/infer_hydrogens.jl")

# Load your CIF file
structure = read("public/outputs/1dht.cif", MMCIFFormat)
model = structure[1]

# Infer hydrogens for backbone N and side chain O
hydrogens = infer_hydrogens(model)

# Print results
println("Inferred ", length(hydrogens), " hydrogens")
for (atom, h_coords) in first(collect(hydrogens), 5)
    println("Atom: ", atomname(atom), " (", resname(atom), ") at ", coords(atom), " -> H at ", h_coords)
end