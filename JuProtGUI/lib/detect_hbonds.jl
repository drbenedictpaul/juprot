# detect_hbonds.jl
using BioStructures
using LinearAlgebra
using Main.JuProtGUI.Utils: hbond_distance


function is_donor(atom)
    atom_name = atomname(atom)
    return startswith(atom_name, "N") && !startswith(atom_name, "NE2")  # Nitrogen donors (exclude histidine NE2)
end

function is_acceptor(atom)
    atom_name = atomname(atom)
    return startswith(atom_name, "O") || (startswith(atom_name, "N") && startswith(atom_name, "NE2"))  # Oxygen or histidine NE2 acceptors
end

function detect_hbonds(protein_atoms, ligand_atoms)
    hbonds = []
    distance_threshold = 3.5

    for p_atom in protein_atoms
        for l_atom in ligand_atoms
            dist = hbond_distance(p_atom, l_atom)
            if dist < distance_threshold
                if is_donor(p_atom) && is_acceptor(l_atom)
                    push!(hbonds, (p_atom, l_atom, dist, "Protein Donor -> Ligand Acceptor"))
                elseif is_donor(l_atom) && is_acceptor(p_atom)
                    push!(hbonds, (p_atom, l_atom, dist, "Ligand Donor -> Protein Acceptor"))
                end
            end
        end
    end
    return hbonds
end

function print_hbonds(hbonds)
    println("\nNumber of hydrogen bonds: ", length(hbonds))
    for hbond in hbonds[1:min(5, length(hbonds))]
        p_atom, l_atom, dist, bond_type = hbond
        println("H-bond: Protein atom: ", atomname(p_atom), " (", resname(p_atom), " ", resnumber(p_atom), ")",
                " - Ligand atom: ", atomname(l_atom), " (", resname(l_atom), ")",
                " - Distance: ", round(dist, digits=2), " Å, Type: ", bond_type)
    end
end