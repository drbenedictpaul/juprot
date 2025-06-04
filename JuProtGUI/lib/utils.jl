module Utils
using LinearAlgebra, BioStructures
export hbond_distance, nonbonded_distance, pipi_distance
function hbond_distance(atom1, atom2)
    return norm(BioStructures.coords(atom1) - BioStructures.coords(atom2))
end
function nonbonded_distance(atom1, atom2)
    return norm(BioStructures.coords(atom1) - BioStructures.coords(atom2))
end
function pipi_distance(atom1, atom2)
    return norm(BioStructures.coords(atom1) - BioStructures.coords(atom2))
end
end