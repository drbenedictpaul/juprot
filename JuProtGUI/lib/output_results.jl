using CSV
   using DataFrames
   using Plots

   function compare_and_save_results(complex_results, output_file, plot_file)
       if length(complex_results) != 2
           println("Error: Exactly two complexes required for comparison.")
           return
       end

       first_name = "first_complex.pdb"
       second_name = "second_complex.pdb"
       if !haskey(complex_results, first_name) || !haskey(complex_results, second_name)
           println("Error: Complex names not found.")
           return
       end
       first_result = complex_results[first_name]
       second_result = complex_results[second_name]

       first_ligand = first_result[:ligand_resname]
       second_ligand = second_result[:ligand_resname]

       first_total = length(first_result[:close_contacts])
       second_total = length(second_result[:close_contacts])
       first_hbonds = length(first_result[:hbonds])
       second_hbonds = length(second_result[:hbonds])

       total_diff = second_total - first_total
       hbond_diff = second_hbonds - first_hbonds

       first_sorted = sort(collect(first_result[:interacting_residues]), by=x->x[2], rev=true)
       second_sorted = sort(collect(second_result[:interacting_residues]), by=x->x[2], rev=true)
       first_top = first_sorted[1:min(3, length(first_sorted))]
       second_top = second_sorted[1:min(3, length(second_sorted))]
       first_top_str = join(["$res ($count)" for (res, count) in first_top], ", ")
       second_top_str = join(["$res ($count)" for (res, count) in second_top], ", ")

       all_residues = union(keys(first_result[:interacting_residues]), keys(second_result[:interacting_residues]))
       changed_residues = String[]
       for res in all_residues
           first_count = get(first_result[:interacting_residues], res, 0)
           second_count = get(second_result[:interacting_residues], res, 0)
           if first_count != second_count
               push!(changed_residues, "$res: First=$first_count, Second=$second_count")
           end
       end
       changed_residues_str = join(changed_residues, "; ")

       first_hbond_res = Set{String}()
       for (p_atom, _, _, _) in first_result[:hbonds]
           push!(first_hbond_res, "$(p_atom.resn) $(p_atom.resi)")
       end
       second_hbond_res = Set{String}()
       for (p_atom, _, _, _) in second_result[:hbonds]
           push!(second_hbond_res, "$(p_atom.resn) $(p_atom.resi)")
       end
       first_hbond_str = join(first_hbond_res, ", ")
       second_hbond_str = join(second_hbond_res, ", ")

       df = DataFrame(
           Metric = [
               "Ligand",
               "Total H-Bonds",
               "H-Bond Difference",
               "Top Interacting Residues",
               "Changed Residues",
               "H-Bond Residues"
           ],
           First_Complex = [
               first_ligand,
               first_hbonds,
               "",
               first_top_str,
               "",
               first_hbond_str
           ],
           Second_Complex = [
               second_ligand,
               second_hbonds,
               hbond_diff,
               second_top_str,
               changed_residues_str,
               second_hbond_str
           ]
       )

       CSV.write(output_file, df, missingstring="")
       println("Side-by-side comparison saved to $output_file")

       detail_rows = []
       for (complex_name, result) in [(first_name, first_result), (second_name, second_result)]
           for (p_atom, _, dist, bond_type) in result[:close_contacts]
               push!(detail_rows, (
                   Complex = complex_name == first_name ? "First Complex" : "Second Complex",
                   Interaction_Type = bond_type,
                   Protein_Residue = "$(p_atom.resn) $(p_atom.resi)",
                   Protein_Atom = p_atom.atom,
                   Ligand_Residue = result[:ligand_resname],
                   Ligand_Atom = p_atom.ligand_atom,
                   Distance = round(dist, digits=2),
                   Angle = round(p_atom.angle, digits=2)
               ))
           end
       end
       detail_df = DataFrame(detail_rows)
       CSV.write(joinpath("public", "outputs", "detailed_interactions.csv"), detail_df)
       println("Detailed interactions saved to public/outputs/detailed_interactions.csv")

       residues = sort(collect(all_residues))
       first_counts = [get(first_result[:interacting_residues], res, 0) for res in residues]
       second_counts = [get(second_result[:interacting_residues], res, 0) for res in residues]

       n = length(residues)
       x = 1:n
       bar_width = 0.35
       plt = bar(
           x .- bar_width/2, first_counts,
           label="First Complex",
           color=:blue,
           alpha=0.6,
           bar_width=bar_width,
           xticks=(1:n, residues),
           xrotation=45,
           xlabel="Residue",
           ylabel="Interaction Count",
           title="Residue Interaction Counts",
           legend=:topright,
           size=(800, 600)
       )
       bar!(plt,
           x .+ bar_width/2, second_counts,
           label="Second Complex",
           color=:red,
           alpha=0.6,
           bar_width=bar_width
       )

       savefig(plt, plot_file)
       println("Bar chart saved to $plot_file")
   end

   function print_analytical_summary(complex_results)
       println("Generating analytical summary for complexes: ", keys(complex_results))
       ordered_complexes = ["first_complex.pdb", "second_complex.pdb"]
       for complex_name in ordered_complexes
           if haskey(complex_results, complex_name)
               result = complex_results[complex_name]
               display_name = complex_name == "first_complex.pdb" ? "First Complex" : "Second Complex"
               println("\nAnalytical Summary for $display_name (Ligand: $(result[:ligand_resname])):")
               println("Total H-Bonds: ", length(result[:close_contacts]))
               sorted_residues = sort(collect(result[:interacting_residues]), by=x->x[2], rev=true)
               top_residues = sorted_residues[1:min(3, length(sorted_residues))]
               println("Top Interacting Residues: ", join(["$res ($count)" for (res, count) in top_residues], ", "))
               hbond_residues = Set{String}()
               for (p_atom, _, _, _) in result[:hbonds]
                   push!(hbond_residues, "$(p_atom.resn) $(p_atom.resi)")
               end
               println("Residues with H-Bonds: ", join(hbond_residues, ", "))
           else
               println("Warning: Complex $complex_name not found in results")
           end
       end
   end