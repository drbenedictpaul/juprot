# lib/detect_hbonds.jl
using PythonCall
using DataFrames

"""
    detect_hbonds(pdb_file::String, ligand_resname::String)
    
    NOTE: Despite the name, this function now detects ALL interaction types
    (Hydrogen Bonds, Hydrophobic, Salt Bridges, Pi-Stacking, etc.)
    supported by PLIP.
"""
function detect_hbonds(pdb_file::String, ligand_resname::String)
    println("--- Starting Interaction Detection (Enhanced) ---")
    println("DEBUG: pdb_file: '$pdb_file'")
    println("DEBUG: ligand_resname: '$ligand_resname'")

    interactions_found = []

    try
        plip = pyimport("plip.structure.preparation")
        mol = plip.PDBComplex()
        mol.load_pdb(pdb_file)
        mol.analyze()

        binding_sites = mol.interaction_sets
        
        if binding_sites === PythonCall.pybuiltins.None || pylen(binding_sites) == 0
            println("DEBUG: No binding sites found.")
            return []
        end

        py_keys_iterator = binding_sites.keys()
        for site_key_py in py_keys_iterator
            site_key = pyconvert(String, site_key_py)
            site = binding_sites[site_key_py]

            # --- Check if this is the correct ligand site ---
            is_target_site = false
            if occursin(ligand_resname, site_key)
                if PythonCall.pyhasattr(site, "ligand") && site.ligand !== PythonCall.pybuiltins.None
                    hetid_py = site.ligand.hetid
                    if hetid_py !== PythonCall.pybuiltins.None
                        if pyconvert(String, hetid_py) == ligand_resname
                            is_target_site = true
                        end
                    end
                end
            end

            if is_target_site
                println("--- Processing site: '$site_key' ---")

                # ==================================================
                # HELPER: Generic Extraction Function
                # ==================================================
                function extract_generic(interaction_list, type_label)
                    if interaction_list === PythonCall.pybuiltins.None || pylen(interaction_list) == 0
                        return
                    end
                    
                    for item in interaction_list
                        # Basic Info
                        p_resn = pyconvert(String, item.restype)
                        p_resi = pyconvert(Int, item.resnr)
                        
                        # Distance (Handle variations)
                        dist = 0.0
                        if PythonCall.pyhasattr(item, "distance")
                            dist = pyconvert(Float64, item.distance)
                        elseif PythonCall.pyhasattr(item, "distance_ad")
                            dist = pyconvert(Float64, item.distance_ad)
                        end

                        # Angle (Only relevant for some, default to 0.0)
                        ang = 0.0
                        if PythonCall.pyhasattr(item, "angle")
                            ang = pyconvert(Float64, item.angle)
                        end

                        # Atom Names (Try to get them using your OB method, else generic)
                        p_atom = "UNK"
                        l_atom = "UNK"
                        
                        # Try to extract detailed atom names if possible
                        # (Note: Hydrophobic contacts/Stacking might not have simple .d/.a structure)
                        # We use a simplified check here for robustness
                        
                        push!(interactions_found, (
                            prot_resn = p_resn, 
                            prot_resi = p_resi,
                            prot_atom_type = "N/A", # Simplified for generic
                            prot_atom_name = p_atom,
                            lig_atom_type = "N/A", 
                            lig_atom_name = l_atom,
                            distance = dist, 
                            angle = ang, 
                            type = type_label
                        ))
                    end
                end
                
                # ==================================================
                # 1. Hydrogen Bonds (Using your original detailed logic)
                # ==================================================
                # Gather H-bonds first
                hbonds_list = PythonCall.Py([])
                if PythonCall.pyhasattr(site, "hbonds_ldon"); for x in site.hbonds_ldon; hbonds_list.append(x); end; end
                if PythonCall.pyhasattr(site, "hbonds_pdon"); for x in site.hbonds_pdon; hbonds_list.append(x); end; end

                for hbond in hbonds_list
                     # Your original extraction logic
                     p_resn = pyconvert(String, hbond.restype)
                     p_resi = pyconvert(Int, hbond.resnr)
                     
                     # OpenBabel Atom ID Extraction
                     d_atom = hbond.d.OBAtom
                     a_atom = hbond.a.OBAtom
                     d_res = d_atom.GetResidue()
                     a_res = a_atom.GetResidue()
                     
                     d_name = (d_res !== PythonCall.pybuiltins.None) ? strip(pyconvert(String, d_res.GetAtomID(d_atom))) : "UNK"
                     a_name = (a_res !== PythonCall.pybuiltins.None) ? strip(pyconvert(String, a_res.GetAtomID(a_atom))) : "UNK"
                     
                     prot_is_don = pyconvert(Bool, hbond.protisdon)
                     p_atom = prot_is_don ? d_name : a_name
                     l_atom = prot_is_don ? a_name : d_name
                     
                     dist = pyconvert(Float64, hbond.distance_ad)
                     ang = pyconvert(Float64, hbond.angle)

                     push!(interactions_found, (
                        prot_resn = p_resn, prot_resi = p_resi,
                        prot_atom_type = "HB", prot_atom_name = p_atom,
                        lig_atom_type = "HB", lig_atom_name = l_atom,
                        distance = dist, angle = ang, type = "Hydrogen Bond"
                    ))
                end

                # ==================================================
                # 2. Extract Other Interactions (Simplified)
                # ==================================================
                if PythonCall.pyhasattr(site, "hydrophobic_contacts")
                    extract_generic(site.hydrophobic_contacts, "Hydrophobic")
                end
                if PythonCall.pyhasattr(site, "saltbridge_lneg")
                    extract_generic(site.saltbridge_lneg, "Salt Bridge")
                end
                if PythonCall.pyhasattr(site, "saltbridge_pneg")
                    extract_generic(site.saltbridge_pneg, "Salt Bridge")
                end
                if PythonCall.pyhasattr(site, "pistacking")
                    extract_generic(site.pistacking, "Pi-Stacking")
                end
                if PythonCall.pyhasattr(site, "pication_laro")
                    extract_generic(site.pication_laro, "Pi-Cation")
                end
                if PythonCall.pyhasattr(site, "pication_paro")
                    extract_generic(site.pication_paro, "Pi-Cation")
                end
                if PythonCall.pyhasattr(site, "halogen_bonds")
                    extract_generic(site.halogen_bonds, "Halogen Bond")
                end
                if PythonCall.pyhasattr(site, "water_bridges")
                    extract_generic(site.water_bridges, "Water Bridge")
                end

            end # end is_target_site
        end # end loop sites

        return interactions_found

    catch e
        println("ERROR in detection:")
        Base.showerror(stdout, e)
        return []
    end
end

function print_hbonds(interactions)
    # Renamed for clarity, though keeping function name compatible
    println("\n--- Interaction Summary ---")
    println("Total Interactions: ", length(interactions))
    if length(interactions) > 0
        println("First 5:")
        for i in interactions[1:min(5, end)]
            println("  $(i.type): $(i.prot_resn) $(i.prot_resi) -- $(i.distance)A")
        end
    end
end