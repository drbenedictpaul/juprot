using BioStructures
using LinearAlgebra
using JuProtGUI.Utils: nonbonded_distance

function is_hydrophobic(atom)
    atom_name = atomname(atom)
    return startswith(atom_name, "C")
end

function detect_nonbonded(protein_atoms, ligand_atoms, min_distance::Float64, max_distance::Float64)
    nonbonded = []
    for p_atom in protein_atoms
        for l_atom in ligand_atoms
            if is_hydrophobic(p_atom) || is_hydrophobic(l_atom)
                dist = nonbonded_distance(p_atom, l_atom)
                if min_distance <= dist <= max_distance
                    push!(nonbonded, (p_atom, l_atom, dist, "Non-Bonded"))
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