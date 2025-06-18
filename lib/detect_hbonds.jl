# lib/detect_hbonds.jl
using PythonCall

"""
    detect_hbonds(pdb_file::String, ligand_resname::String)
... (docstring) ...
"""
function detect_hbonds(pdb_file::String, ligand_resname::String)
    println("--- Starting H-bond detection in detect_hbonds ---")
    println("DEBUG: pdb_file: '$pdb_file'")
    println("DEBUG: ligand_resname passed to detect_hbonds is: '$ligand_resname'")

    hbonds_found = []
    try
        plip = pyimport("plip.structure.preparation")
        # ob_py = pyimport("openbabel") # Not strictly needed to import if just calling methods

        mol = plip.PDBComplex()
        mol.load_pdb(pdb_file)
        println("DEBUG: PDB file loaded.")
        mol.analyze()
        println("DEBUG: mol.analyze() called.")

        binding_sites = mol.interaction_sets
        
        if binding_sites === PythonCall.pybuiltins.None || pylen(binding_sites) == 0
            println("DEBUG: No binding sites found by PLIP.")
            return []
        else
            println("DEBUG: PLIP found $(pylen(binding_sites)) binding sites.")
        end

        py_keys_iterator = binding_sites.keys()
        for site_key_py in py_keys_iterator
            site_key = pyconvert(String, site_key_py)
            site = binding_sites[site_key_py]
            println("DEBUG: Checking site_key: '$site_key'")

            if occursin(ligand_resname, site_key)
                if PythonCall.pyhasattr(site, "ligand") && site.ligand !== PythonCall.pybuiltins.None && PythonCall.pyhasattr(site.ligand, "hetid")
                    plip_ligand_hetid_py = site.ligand.hetid
                    if plip_ligand_hetid_py !== PythonCall.pybuiltins.None
                        plip_ligand_hetid = pyconvert(String, plip_ligand_hetid_py)
                        if plip_ligand_hetid == ligand_resname
                            println("--- Processing site: '$site_key' for target ligand: '$ligand_resname' (Exact HETID Match) ---")
                            all_site_hbonds_py = PythonCall.Py([])
                            
                            if PythonCall.pyhasattr(site, "hbonds_ldon") && site.hbonds_ldon !== PythonCall.pybuiltins.None && pylen(site.hbonds_ldon) > 0
                                println("DEBUG: Found $(pylen(site.hbonds_ldon)) H-bonds in site.hbonds_ldon.")
                                for item in site.hbonds_ldon; all_site_hbonds_py.append(item); end
                            else
                                println("DEBUG: No/empty site.hbonds_ldon.")
                            end
                            if PythonCall.pyhasattr(site, "hbonds_pdon") && site.hbonds_pdon !== PythonCall.pybuiltins.None && pylen(site.hbonds_pdon) > 0
                                println("DEBUG: Found $(pylen(site.hbonds_pdon)) H-bonds in site.hbonds_pdon.")
                                for item in site.hbonds_pdon; all_site_hbonds_py.append(item); end
                            else
                                 println("DEBUG: No/empty site.hbonds_pdon.")
                            end

                            println("DEBUG: Total $(pylen(all_site_hbonds_py)) H-bonds to process for site '$site_key'.")
                            
                            for hbond_py in all_site_hbonds_py
                                println("DEBUG Hbond object from PLIP: $(pystr(hbond_py))")

                                local protein_resn, protein_resi, protein_atom_name, protein_atom_type
                                local ligand_atom_name, ligand_atom_type

                                protein_resn = pyconvert(String, hbond_py.restype)
                                protein_resi = pyconvert(Int, hbond_py.resnr)

                                pyb_d_atom = hbond_py.d 
                                pyb_a_atom = hbond_py.a 

                                ob_d_atom = pyb_d_atom.OBAtom
                                ob_a_atom = pyb_a_atom.OBAtom

                                ob_d_residue = ob_d_atom.GetResidue()
                                ob_a_residue = ob_a_atom.GetResidue()
                                
                                d_atom_pdb_name_str = "UNK_D"
                                a_atom_pdb_name_str = "UNK_A"

                                if ob_d_residue !== PythonCall.pybuiltins.None
                                    d_atom_pdb_name_py = ob_d_residue.GetAtomID(ob_d_atom)
                                    d_atom_pdb_name_str = strip(pyconvert(String, d_atom_pdb_name_py))
                                else
                                    println("WARN: Donor atom for hbond $(pystr(hbond_py)) has no OBResidue.")
                                end
                                if ob_a_residue !== PythonCall.pybuiltins.None
                                    a_atom_pdb_name_py = ob_a_residue.GetAtomID(ob_a_atom)
                                    a_atom_pdb_name_str = strip(pyconvert(String, a_atom_pdb_name_py))
                                else
                                     println("WARN: Acceptor atom for hbond $(pystr(hbond_py)) has no OBResidue.")
                                end
                                
                                println("DEBUG: Donor PDB Name: '$(d_atom_pdb_name_str)', Acceptor PDB Name: '$(a_atom_pdb_name_str)'")

                                if pyconvert(Bool, hbond_py.protisdon) 
                                    protein_atom_name = d_atom_pdb_name_str
                                    protein_atom_type = pyconvert(String, hbond_py.dtype) 
                                    
                                    ligand_atom_name = a_atom_pdb_name_str
                                    ligand_atom_type = pyconvert(String, hbond_py.atype)  
                                else 
                                    protein_atom_name = a_atom_pdb_name_str
                                    protein_atom_type = pyconvert(String, hbond_py.atype)
                                    
                                    ligand_atom_name = d_atom_pdb_name_str
                                    ligand_atom_type = pyconvert(String, hbond_py.dtype)
                                end
                                
                                distance = pyconvert(Float64, hbond_py.distance_ad)
                                angle = pyconvert(Float64, hbond_py.angle)

                                interaction = (
                                    prot_resn = protein_resn, prot_resi = protein_resi,
                                    prot_atom_type = protein_atom_type, prot_atom_name = protein_atom_name,
                                    lig_atom_type = ligand_atom_type, lig_atom_name = ligand_atom_name,
                                    distance = distance, angle = angle, type = "Hydrogen Bond"
                                )
                                push!(hbonds_found, interaction)
                                println("DEBUG: Protein-Ligand H-bond ADDED: Prot $(protein_atom_name)($(protein_resn) $(protein_resi)) - Lig $(ligand_atom_name)($(ligand_resname)), D: $(round(distance, digits=2)) Å, A: $(round(angle, digits=2))°")
                            end 
                        end # end if HETID match
                    end # end if hetid_py not None
                end # end if hasattr ligand and hetid
            end # end if occursin
            println("--- Finished checking site: '$site_key' ---")
        end # End loop over binding_sites

        println("Total H-bonds found for ligand '$ligand_resname' in '$pdb_file' after all filtering: ", length(hbonds_found))
        return hbonds_found
    catch e 
        println("ERROR: An error occurred in detect_hbonds for '$pdb_file', ligand '$ligand_resname':")
        Base.showerror(stdout, e) 
        Base.show_backtrace(stdout, Base.catch_backtrace()) 
        println() 

        if e isa PythonCall.PyException
            println("--- Wrapped Python Exception Details ---")
            py_exc_obj = e.exc 
            println("Python Exception (pyrepr): $(PythonCall.pyrepr(String, py_exc_obj))")
            py_traceback_mod = pyimport("traceback")
            println("Attempting to print Python traceback using traceback.print_exception():")
            try
                py_exc_type = pytype(py_exc_obj)
                py_exc_tb = PythonCall.pyhasattr(py_exc_obj, "__traceback__") ? py_exc_obj.__traceback__ : PythonCall.pybuiltins.None
                py_traceback_mod.print_exception(py_exc_type, py_exc_obj, py_exc_tb)
            catch py_print_err
                println("Failed to print Python traceback using traceback.print_exception(): $py_print_err")
                try # Fallback
                    formatted_tb = py_traceback_mod.format_exc()
                    println("Fallback format_exc():\n", pyconvert(String, formatted_tb))
                catch forma_exc_err
                    println("format_exc() also failed: $forma_exc_err")
                end
            end
            println("------------------------------------")
        else
            println("The caught exception was a Julia exception, not directly a Python one.")
        end
        return [] 
    end
end

function print_hbonds(hbonds)
    count = length(hbonds)
    println("\n--- Hydrogen Bond Summary (from print_hbonds function) ---")
    println("Number of hydrogen bonds: ", count)
    if count > 0
        println("\nShowing details for up to the first 5 bonds:")
        for bond in hbonds[1:min(5, count)]
            println("  - Protein: ", bond.prot_atom_name, " (", bond.prot_resn, " ", bond.prot_resi, ")",
                    " -- Ligand: ", bond.lig_atom_name,
                    " | Distance: ", round(bond.distance, digits=2), " Å",
                    " | Angle: ", round(bond.angle, digits=2), "°")
        end
    end
    println("---------------------------\n")
end