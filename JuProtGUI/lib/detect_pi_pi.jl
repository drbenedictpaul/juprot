# detect_pi_pi.jl
using BioStructures
using LinearAlgebra
using JuProtGUI.Utils: pipi_distance

function is_aromatic(atom)
    res_name = resname(atom)
    return res_name in ["PHE", "TYR", "TRP", "HIS"] || res_name == "ASD"  # Aromatic residues or ASD
end

function detect_pi_pi(protein_atoms, ligand_atoms)
    pi_pi = []
    distance_threshold = 6.0  # Å for π-π stacking

    for p_atom in protein_atoms
        if is_aromatic(p_atom)
            for l_atom in ligand_atoms
                if is_aromatic(l_atom)
                    dist = pipi_distance(p_atom, l_atom)
                    if dist < distance_threshold
                        push!(pi_pi, (p_atom, l_atom, dist, "π-π Stacking"))
                    end
                end
            end
        end
    end
    return pi_pi
end

function print_pi_pi(pi_pi)
    println("\nNumber of π-π stacking interactions: ", length(pi_pi))
    for interaction in pi_pi[1:min(5, length(pi_pi))]
        p_atom, l_atom, dist, bond_type = interaction
        println("π-π: Protein atom: ", atomname(p_atom), " (", resname(p_atom), " ", resnumber(p_atom), ")",
                " - Ligand atom: ", atomname(l_atom), " (", resname(l_atom), ")",
                " - Distance: ", round(dist, digits=2), " Å, Type: ", bond_type)
    end
end
