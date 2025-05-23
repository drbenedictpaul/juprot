# parse_pdb.jl
using BioStructures
using LinearAlgebra

# Include detection and output scripts
include("detect_hbonds.jl")
include("detect_nonbonded.jl")
include("output_results.jl")

# Prompt user for exactly two CIF files
println("Enter the native CIF file path:")
native_cif = strip(readline())
println("Enter the mutated CIF file path:")
mutated_cif = strip(readline())
cif_files = [native_cif, mutated_cif]

if any(isempty, cif_files)
    println("Error: Both native and mutated CIF files must be provided. Exiting.")
    exit(1)
end

# Store results and structures for each complex
complex_results = Dict{String, Dict}()
structures = Dict{String, Any}()

# Process each CIF file
for cif_file in cif_files
    println("\nProcessing $cif_file...")
    # Read the CIF file
    structure = read(cif_file, MMCIFFormat)
    models = collectmodels(structure)
    println("Number of models: ", length(models))
    model = structure[1]

    # Select protein atoms
    protein_atoms = collectatoms(model, standardselector)
    println("Number of protein atoms: ", length(protein_atoms))

    # List non-protein, non-water residue names
    non_protein_atoms = collectatoms(model, atom -> !standardselector(atom) && resname(atom) != "HOH" && resname(atom) != "WAT")
    unique_resnames = unique(resname(atom) for atom in non_protein_atoms)
    println("Available ligand residue names: ", unique_resnames)

    # Prompt user to select a ligand
    println("Enter the ligand residue name to analyze for $cif_file (e.g., ASD for androstenedione): ")
    ligand_resname = strip(readline())
    if !(ligand_resname in unique_resnames)
        println("Error: '$ligand_resname' not found in available ligands. Skipping $cif_file.")
        continue
    end

    # Define ligand selector
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

    # Identify interacting amino acids
    interacting_residues = Dict{String, Int}()
    for (p_atom, l_atom, dist) in close_contacts
        res_key = "$(resname(p_atom)) $(resnumber(p_atom))"
        interacting_residues[res_key] = get(interacting_residues, res_key, 0) + 1
    end

    # Find residue with maximum interactions
    max_residue = ""
    max_count = 0
    for (res, count) in interacting_residues
        if count > max_count
            max_residue = res
            max_count = count
        end
    end

    # Detect hydrogen bonds and non-bonded interactions
    hbonds = detect_hbonds(protein_atoms, ligand_atoms)
    nonbonded = detect_nonbonded(protein_atoms, ligand_atoms)

    # Store results and structure
    complex_name = basename(cif_file)
    complex_results[complex_name] = Dict(
        :close_contacts => close_contacts,
        :hbonds => hbonds,
        :nonbonded => nonbonded,
        :interacting_residues => interacting_residues,
        :max_residue => max_residue,
        :max_count => max_count,
        :ligand_resname => ligand_resname
    )
    structures[complex_name] = structure
end

# Generate and save side-by-side comparison
compare_and_save_results(complex_results, "comparison_table.csv", "residue_interactions.png", 
                        structures[cif_files[1]], structures[cif_files[2]])
print_analytical_summary(complex_results)