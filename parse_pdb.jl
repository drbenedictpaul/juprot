# Import the required packages
using PDBTools

# Specify the path to your PDB file
pdb_file = "3eqm.cif"  # Replace with the actual path to your PDB file

# Read the PDB file
atoms = readPDB(pdb_file)

# Print the total number of atoms to verify loading
println("Total number of atoms: ", length(atoms))

# Define a selection for protein atoms (standard amino acid residues)
protein = select(atoms, "protein")
println("Number of protein atoms: ", length(protein))

# Define a selection for ligand atoms (not protein, not water)
ligand = select(atoms, "not protein and not water")
println("Number of ligand atoms: ", length(ligand))

# Optional: Print some details about the ligand atoms
for atom in ligand[1:min(5, length(ligand))]  # Print first 5 ligand atoms
    println("Atom: ", atom.name, " Residue: ", atom.resname, " Chain: ", atom.chain)
end