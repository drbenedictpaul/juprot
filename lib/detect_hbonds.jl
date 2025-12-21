using PythonCall
using DataFrames

function detect_hbonds(pdb_file::String, ligand_resname::String)
    interactions_found = []
    ligand_meta = Dict("chain" => "", "resi" => 0)
    
    try
        # Standard PLIP import
        plip_prep = pyimport("plip.structure.preparation")
        
        mol = plip_prep.PDBComplex()
        mol.load_pdb(pdb_file)
        mol.analyze()
        
        binding_sites = mol.interaction_sets
        
        if binding_sites === PythonCall.pybuiltins.None || pylen(binding_sites) == 0
            return (interactions_found, ligand_meta)
        end

        for (site_key_py, site) in pyconvert(Dict, binding_sites)
            site_key = pyconvert(String, site_key_py)
            
            if occursin(ligand_resname, site_key) && PythonCall.pyhasattr(site, "ligand") && 
               site.ligand !== PythonCall.pybuiltins.None && pyconvert(String, site.ligand.hetid) == ligand_resname
                
                # Capture Ligand Location
                lig_chain = PythonCall.pyhasattr(site.ligand, "chain") ? pyconvert(String, site.ligand.chain) : ""
                lig_resi  = PythonCall.pyhasattr(site.ligand, "resnr") ? pyconvert(Int, site.ligand.resnr) : 0
                ligand_meta = Dict("chain" => lig_chain, "resi" => lig_resi)

                # Extract Interactions
                interaction_sources = [
                    (PythonCall.pyhasattr(site, "hbonds_pdon") ? site.hbonds_pdon : [], "Hydrogen Bond"),
                    (PythonCall.pyhasattr(site, "hbonds_ldon") ? site.hbonds_ldon : [], "Hydrogen Bond"),
                    (PythonCall.pyhasattr(site, "hydrophobic_contacts") ? site.hydrophobic_contacts : [], "Hydrophobic"),
                    (PythonCall.pyhasattr(site, "saltbridge_lneg") ? site.saltbridge_lneg : [], "Salt Bridge"),
                    (PythonCall.pyhasattr(site, "saltbridge_pneg") ? site.saltbridge_pneg : [], "Salt Bridge"),
                    (PythonCall.pyhasattr(site, "pistacking") ? site.pistacking : [], "Pi-Stacking"),
                    (PythonCall.pyhasattr(site, "pication_laro") ? site.pication_laro : [], "Pi-Cation"),
                    (PythonCall.pyhasattr(site, "pication_paro") ? site.pication_paro : [], "Pi-Cation"),
                    (PythonCall.pyhasattr(site, "halogen_bonds") ? site.halogen_bonds : [], "Halogen Bond"),
                    (PythonCall.pyhasattr(site, "water_bridges") ? site.water_bridges : [], "Water Bridge")
                ]

                for (interaction_list, type_label) in interaction_sources
                    if !isempty(interaction_list) && interaction_list !== PythonCall.pybuiltins.None
                        for item in interaction_list
                            p_chain = PythonCall.pyhasattr(item, "reschain") ? pyconvert(String, item.reschain) : ""
                            push!(interactions_found, (
                                prot_resn = pyconvert(String, item.restype), 
                                prot_resi = pyconvert(Int, item.resnr), 
                                prot_chain = p_chain,
                                type = type_label
                            ))
                        end
                    end
                end
                break 
            end
        end
    catch e
        println("ERROR in interaction detection: ")
        Base.showerror(stdout, e)
    end
    
    return (interactions_found, ligand_meta)
end