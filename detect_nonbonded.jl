# detect_nonbonded.jl
using BioStructures
using LinearAlgebra

function distance(atom1, atom2)
    return norm(BioStructures.coords(atom1) - BioStructures.coords(atom2))
end

function is_hydrophobic(atom)
    atom_name = atomname(atom)
    return startswith(atom_name, "C")
end

function detect_nonbonded(protein_atoms, ligand_atoms)
    nonbonded = []
    distance_threshold = 4.0  # Align with close contact threshold

    for p_atom in protein_atoms
        if is_hydrophobic(p_atom)
            for l_atom in ligand_atoms
                if is_hydrophobic(l_atom)
                    dist = distance(p_atom, l_atom)
                    if dist < distance_threshold
                        push!(nonbonded, (p_atom, l_atom, dist, "Hydrophobic Interaction"))
                    end
                end
            end
        end
    end
    return nonbonded
end

function print_nonbonded(nonbonded)
    println("\nNumber of non-bonded (hydrophobic) interactions: ", length(nonbonded))
    for nb in nonbonded[1:min(5, length(nonbonded))]
        p_atom, l_atom, dist, bond_type = nb
        println("Non-bonded: Protein atom: ", atomname(p_atom), " (", resname(p_atom), " ", resnumber(p_atom), ")",
                " - Ligand atom: ", atomname(l_atom), " (", resname(l_atom), ")",
                " - Distance: ", round(dist, digits=2), " Å, Type: ", bond_type)
    end
end