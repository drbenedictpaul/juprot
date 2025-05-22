# Import necessary packages
using BioStructures
using LinearAlgebra  # For distance calculations

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

# Define a custom ligand selector (non-protein, non-water, e.g., ACT for 3EQM)
function ligand_selector(atom)
    !standardselector(atom) && resname(atom) != "HOH" && resname(atom) != "WAT" && resname(atom) == "ACT"
end
ligand_atoms = collectatoms(model, ligand_selector)
println("Number of ligand atoms: ", length(ligand_atoms))

# Function to calculate distance between two atoms
function distance(atom1, atom2)
    return norm(coords(atom1) - coords(atom2))
end

# Set a distance threshold for interactions (e.g., 4.0 Å)
threshold = 4.0

# Find close contacts between protein and ligand atoms
close_contacts = []
for p_atom in protein_atoms
    for l_atom in ligand_atoms
        dist = distance(p_atom, l_atom)
        if dist < threshold
            push!(close_contacts, (p_atom, l_atom, dist))
        end
    end
end

# Print the number of close contacts found
println("Number of close contacts: ", length(close_contacts))

# Print details of the first few close contacts
for contact in close_contacts[1:min(5, length(close_contacts))]
    p_atom, l_atom, dist = contact
    println("Protein atom: ", atomname(p_atom), " (", resname(p_atom), " ", resnumber(p_atom), ")",
            " - Ligand atom: ", atomname(l_atom), " (", resname(l_atom), ")",
            " - Distance: ", round(dist, digits=2), " Å")
end

# Optional: Inspect ligand residue names to verify selection
if !isempty(ligand_atoms)
    unique_resnames = unique(resname(atom) for atom in ligand_atoms)
    println("Ligand residue names: ", unique_resnames)
else
    println("Warning: No ligand atoms detected.")
end