# parse_pdb.jl
using BioStructures
using LinearAlgebra

# Include the hydrogen bond detection script
include("detect_hbonds.jl")

# Specify the path to your CIF file
cif_file = "3eqm.cif"  # Update with the full path if not in the project directory

# Read the CIF file
structure = read(cif_file, MMCIFFormat)

# Check the number of models
models = collectmodels(structure)
println("Number of models: ", length(models))

# Select the first model
model = structure[1]

# Select protein atoms (standard amino acid residues)
protein_atoms = collectatoms(model, standardselector)
println("Number of protein atoms: ", length(protein_atoms))

# List all non-protein, non-water residue names
non_protein_atoms = collectatoms(model, atom -> !standardselector(atom) && resname(atom) != "HOH" && resname(atom) != "WAT")
unique_resnames = unique(resname(atom) for atom in non_protein_atoms)
println("Available ligand residue names: ", unique_resnames)

# Prompt user to select a ligand
println("Enter the ligand residue name to analyze (e.g., ASD for androstenedione): ")
ligand_resname = strip(readline())
if !(ligand_resname in unique_resnames)
    println("Error: '$ligand_resname' not found in available ligands. Exiting.")
    exit(1)
end

# Define a custom ligand selector based on user input
function ligand_selector(atom)
    resname(atom) == ligand_resname
end
ligand_atoms = collectatoms(model, ligand_selector)
println("Number of ligand atoms: ", length(ligand_atoms))

# Find close contacts
close_contacts = []
threshold = 4.0
for p_atom in protein_atoms
    for l_atom in ligand_atoms
        dist = distance(p_atom, l_atom)
        if dist < threshold
            push!(close_contacts, (p_atom, l_atom, dist))
        end
    end
end

# Identify interacting amino acids and count interactions
interacting_residues = Dict{String, Int}()
for (p_atom, l_atom, dist) in close_contacts
    res_key = "$(resname(p_atom)) $(resnumber(p_atom))"
    interacting_residues[res_key] = get(interacting_residues, res_key, 0) + 1
end

# Find residue with maximum interactions
let
    max_residue = ""
    max_count = 0
    for (res, count) in interacting_residues
        if count > max_count
            max_residue = res
            max_count = count
        end
    end

    # Print close contacts
    println("Number of close contacts: ", length(close_contacts))
    for contact in close_contacts[1:min(5, length(close_contacts))]
        p_atom, l_atom, dist = contact
        println("Close contact: Protein atom: ", atomname(p_atom), " (", resname(p_atom), " ", resnumber(p_atom), ")",
                " - Ligand atom: ", atomname(l_atom), " (", resname(l_atom), ")",
                " - Distance: ", round(dist, digits=2), " Å")
    end

    # Detect and print hydrogen bonds
    hbonds = detect_hbonds(protein_atoms, ligand_atoms)
    print_hbonds(hbonds)

    # Print interacting residues
    println("\nInteracting amino acids and interaction counts:")
    for (res, count) in sort(collect(interacting_residues), by=x->x[2], rev=true)
        println("Residue: $res, Interactions: $count")
    end
    if !isempty(max_residue)
        println("Residue with maximum interactions: $max_residue ($max_count interactions)")
    else
        println("No interacting residues found.")
    end
end

# Inspect ligand residue names
if !isempty(ligand_atoms)
    unique_ligand_resnames = unique(resname(atom) for atom in ligand_atoms)
    println("Selected ligand residue name: ", unique_ligand_resnames)
else
    println("Warning: No ligand atoms detected for '$ligand_resname'.")
end