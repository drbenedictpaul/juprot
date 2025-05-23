# detect_hbonds.jl
using BioStructures
using LinearAlgebra

function distance(atom1, atom2)
    return norm(BioStructures.coords(atom1) - BioStructures.coords(atom2))
end

function is_donor_or_acceptor(atom)
    atom_name = atomname(atom)
    return startswith(atom_name, "N") || startswith(atom_name, "O")
end

function detect_hbonds(protein_atoms, ligand_atoms)
    hbonds = []
    distance_threshold = 3.5

    for p_atom in protein_atoms
        if is_donor_or_acceptor(p_atom)
            for l_atom in ligand_atoms
                if is_donor_or_acceptor(l_atom)
                    dist = distance(p_atom, l_atom)
                    if dist < distance_threshold
                        push!(hbonds, (p_atom, l_atom, dist))
                    end
                end
            end
        end
    end
    return hbonds
end

function print_hbonds(hbonds)
    println("\nNumber of hydrogen bonds: ", length(hbonds))
    for hbond in hbonds[1:min(5, length(hbonds))]
        p_atom, l_atom, dist = hbond
        println("H-bond: Protein atom: ", atomname(p_atom), " (", resname(p_atom), " ", resnumber(p_atom), ")",
                " - Ligand atom: ", atomname(l_atom), " (", resname(l_atom), ")",
                " - Distance: ", round(dist, digits=2), " Å")
    end
end