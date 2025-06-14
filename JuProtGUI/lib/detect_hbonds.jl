using PythonCall

   function detect_hbonds(pdb_file, ligand_resname)
       println("Starting H-bond detection with Python PLIP for ", pdb_file)
       try
           plip = pyimport("plip.structure.preparation")
           mol = plip.PDBComplex()
           mol.load_pdb(pdb_file)
           mol.analyze()
           binding_sites = mol.interaction_sets
           hbonds = []
           for site_key in keys(binding_sites)
               site = binding_sites[site_key]
               if occursin(ligand_resname, site_key)
                   for hbond in site.hbonds
                       donor_res = hbond.d_res
                       donor_num = hbond.d_resno
                       donor_atom = hbond.d_atom
                       acceptor_atom = hbond.a_atom
                       distance = hbond.distance_ad
                       angle = hbond.angle
                       interaction = (
                           resn=donor_res,
                           resi=donor_num,
                           atom=donor_atom,
                           ligand_atom=acceptor_atom,
                           distance=distance,
                           angle=angle
                       )
                       push!(hbonds, (interaction, py"None", distance, "Hydrogen Bond"))
                       println("H-bond detected: Protein ", donor_atom, " (", donor_res, " ", donor_num, ") - Ligand ", acceptor_atom, " (", ligand_resname, "), Distance: ", round(distance, digits=2), " Å, Angle: ", round(angle, digits=2), "°")
                   end
               end
           end
           println("Total filtered H-bonds: ", length(hbonds))
           return hbonds
       catch e
           println("Error running PLIP: ", e)
           return []
       end
   end

   function print_hbonds(hbonds)
       println("\nNumber of hydrogen bonds: ", length(hbonds))
       for hbond in hbonds[1:min(5, length(hbonds))]
           p_atom, _, dist, bond_type = hbond
           println("H-bond: Protein atom: ", p_atom.atom, " (", p_atom.resn, " ", p_atom.resi, ")",
                   " - Ligand atom: ", p_atom.ligand_atom,
                   " - Distance: ", round(dist, digits=2), " Å, Angle: ", round(p_atom.angle, digits=2), "°, Type: ", bond_type)
       end
   end