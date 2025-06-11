# output_results.jl
using CSV
using DataFrames
using Plots
using BioStructures

function check_residue_differences(native_structure, mutated_structure)
    native_model = native_structure[1]
    mutated_model = mutated_structure[1]
    native_residues = Dict{String, String}()
    mutated_residues = Dict{String, String}()

    for res in collectresidues(native_model, standardselector)
        key = "$(resnumber(res))"
        native_residues[key] = resname(res)
    end
    for res in collectresidues(mutated_model, standardselector)
        key = "$(resnumber(res))"
        mutated_residues[key] = resname(res)
    end

    differences = String[]
    for key in union(keys(native_residues), keys(mutated_residues))
        native_res = get(native_residues, key, "None")
        mutated_res = get(mutated_residues, key, "None")
        if native_res != mutated_res
            push!(differences, "Position $key: Native=$native_res, Mutated=$mutated_res")
        end
    end
    return differences
end

function compare_and_save_results(complex_results, output_file, plot_file, native_structure, mutated_structure)
    if length(complex_results) != 2
        println("Error: Exactly two complexes required for comparison.")
        return
    end

    residue_diffs = check_residue_differences(native_structure, mutated_structure)
    if isempty(residue_diffs)
        println("No residue differences detected between native and mutated proteins.")
    else
        println("Residue differences between native and mutated proteins:")
        for diff in residue_diffs
            println(diff)
        end
    end

    first_name = "first_complex.cif"
    second_name = "second_complex.cif"
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
    first_nonbonded = length(first_result[:nonbonded])
    second_nonbonded = length(second_result[:nonbonded])
    first_pi_pi = length(first_result[:pi_pi])
    second_pi_pi = length(second_result[:pi_pi])

    total_diff = second_total - first_total
    hbond_diff = second_hbonds - first_hbonds
    nonbonded_diff = second_nonbonded - first_nonbonded
    pi_pi_diff = second_pi_pi - first_pi_pi

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
        push!(first_hbond_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    second_hbond_res = Set{String}()
    for (p_atom, _, _, _) in second_result[:hbonds]
        push!(second_hbond_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    first_hbond_str = join(first_hbond_res, ", ")
    second_hbond_str = join(second_hbond_res, ", ")

    first_pi_pi_res = Set{String}()
    for (p_atom, _, _, _) in first_result[:pi_pi]
        push!(first_pi_pi_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    second_pi_pi_res = Set{String}()
    for (p_atom, _, _, _) in second_result[:pi_pi]
        push!(second_pi_pi_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    first_pi_pi_str = join(first_pi_pi_res, ", ")
    second_pi_pi_str = join(second_pi_pi_res, ", ")

    first_hbond_pct = first_total > 0 ? round((first_hbonds / first_total) * 100, digits=1) : 0.0
    second_hbond_pct = second_total > 0 ? round((second_hbonds / second_total) * 100, digits=1) : 0.0
    first_nonbonded_pct = first_total > 0 ? round((first_nonbonded / first_total) * 100, digits=1) : 0.0
    second_nonbonded_pct = second_total > 0 ? round((second_nonbonded / second_total) * 100, digits=1) : 0.0
    first_insight = first_hbond_pct > first_nonbonded_pct ? "H-bonds dominate ($first_hbond_pct%)" : "Hydrophobic dominate ($first_nonbonded_pct%)"
    second_insight = second_hbond_pct > second_nonbonded_pct ? "H-bonds dominate ($second_hbond_pct%)" : "Hydrophobic dominate ($second_nonbonded_pct%)"

    df = DataFrame(
        Metric = [
            "Ligand",
            "Total Interactions",
            "Interaction Difference",
            "Hydrogen Bonds",
            "H-Bond Difference",
            "Non-Bonded (Hydrophobic)",
            "Non-Bonded Difference",
            "π-π Stacking",
            "π-π Difference",
            "Top Interacting Residues",
            "Changed Residues",
            "H-Bond Residues",
            "π-π Residues",
            "Binding Insight"
        ],
        First_Complex = [
            first_ligand,
            first_total,
            "",
            first_hbonds,
            "",
            first_nonbonded,
            "",
            first_pi_pi,
            "",
            first_top_str,
            "",
            first_hbond_str,
            first_pi_pi_str,
            first_insight
        ],
        Second_Complex = [
            second_ligand,
            second_total,
            total_diff,
            second_hbonds,
            hbond_diff,
            second_nonbonded,
            nonbonded_diff,
            second_pi_pi,
            pi_pi_diff,
            second_top_str,
            changed_residues_str,
            second_hbond_str,
            second_pi_pi_str,
            second_insight
        ]
    )

    CSV.write(output_file, df, missingstring="")
    println("Side-by-side comparison saved to $output_file")

    detail_rows = []
    for (complex_name, result) in [(first_name, first_result), (second_name, second_result)]
        for (p_atom, l_atom, dist, bond_type) in result[:close_contacts]
            push!(detail_rows, (
                Complex = complex_name == first_name ? "First Complex" : "Second Complex",
                Interaction_Type = "Close Contact",
                Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
                Protein_Atom = atomname(p_atom),
                Ligand_Residue = resname(l_atom),
                Ligand_Atom = atomname(l_atom),
                Distance = round(dist, digits=2),
                Bond_Description = bond_type
            ))
        end
        for (p_atom, l_atom, dist, bond_type) in result[:hbonds]
            push!(detail_rows, (
                Complex = complex_name == first_name ? "First Complex" : "Second Complex",
                Interaction_Type = "Hydrogen Bond",
                Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
                Protein_Atom = atomname(p_atom),
                Ligand_Residue = resname(l_atom),
                Ligand_Atom = atomname(l_atom),
                Distance = round(dist, digits=2),
                Bond_Description = bond_type
            ))
        end
        for (p_atom, l_atom, dist, bond_type) in result[:nonbonded]
            push!(detail_rows, (
                Complex = complex_name == first_name ? "First Complex" : "Second Complex",
                Interaction_Type = "Non-Bonded",
                Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
                Protein_Atom = atomname(p_atom),
                Ligand_Residue = resname(l_atom),
                Ligand_Atom = atomname(l_atom),
                Distance = round(dist, digits=2),
                Bond_Description = bond_type
            ))
        end
        for (p_atom, l_atom, dist, bond_type) in result[:pi_pi]
            push!(detail_rows, (
                Complex = complex_name == first_name ? "First Complex" : "Second Complex",
                Interaction_Type = "π-π Stacking",
                Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
                Protein_Atom = atomname(p_atom),
                Ligand_Residue = resname(l_atom),
                Ligand_Atom = atomname(l_atom),
                Distance = round(dist, digits=2),
                Bond_Description = bond_type
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
        label="First Complex (first_complex.cif)",
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
        label="Second Complex (second_complex.cif)",
        color=:red,
        alpha=0.6,
        bar_width=bar_width
    )

    savefig(plt, plot_file)
    println("Bar chart saved to $plot_file")
end
# function compare_and_save_results(complex_results, output_file, plot_file, native_structure, mutated_structure)
#     if length(complex_results) != 2
#         println("Error: Exactly two complexes (native and mutated) required for comparison.")
#         return
#     end

#     residue_diffs = check_residue_differences(native_structure, mutated_structure)
#     if isempty(residue_diffs)
#         println("No residue differences detected between native and mutated proteins.")
#     else
#         println("Residue differences between native and mutated proteins:")
#         for diff in residue_diffs
#             println(diff)
#         end
#     end

#     complexes = collect(keys(complex_results))
#     native_name, mutated_name = complexes[1], complexes[2]
#     native_result = complex_results[native_name]
#     mutated_result = complex_results[mutated_name]

#     native_ligand = native_result[:ligand_resname]
#     mutated_ligand = mutated_result[:ligand_resname]

#     native_total = length(native_result[:close_contacts])
#     mutated_total = length(mutated_result[:close_contacts])
#     native_hbonds = length(native_result[:hbonds])
#     mutated_hbonds = length(mutated_result[:hbonds])
#     native_nonbonded = length(native_result[:nonbonded])
#     mutated_nonbonded = length(mutated_result[:nonbonded])
#     native_pi_pi = length(native_result[:pi_pi])
#     mutated_pi_pi = length(mutated_result[:pi_pi])

#     total_diff = mutated_total - native_total
#     hbond_diff = mutated_hbonds - native_hbonds
#     nonbonded_diff = mutated_nonbonded - native_nonbonded
#     pi_pi_diff = mutated_pi_pi - native_pi_pi

#     native_sorted = sort(collect(native_result[:interacting_residues]), by=x->x[2], rev=true)
#     mutated_sorted = sort(collect(mutated_result[:interacting_residues]), by=x->x[2], rev=true)
#     native_top = native_sorted[1:min(3, length(native_sorted))]
#     mutated_top = mutated_sorted[1:min(3, length(mutated_sorted))]
#     native_top_str = join(["$res ($count)" for (res, count) in native_top], ", ")
#     mutated_top_str = join(["$res ($count)" for (res, count) in mutated_top], ", ")

#     all_residues = union(keys(native_result[:interacting_residues]), keys(mutated_result[:interacting_residues]))
#     changed_residues = String[]
#     for res in all_residues
#         native_count = get(native_result[:interacting_residues], res, 0)
#         mutated_count = get(mutated_result[:interacting_residues], res, 0)
#         if native_count != mutated_count
#             push!(changed_residues, "$res: Native=$native_count, Mutated=$mutated_count")
#         end
#     end
#     changed_residues_str = join(changed_residues, "; ")

#     native_hbond_res = Set{String}()
#     for (p_atom, _, _, _) in native_result[:hbonds]
#         push!(native_hbond_res, "$(resname(p_atom)) $(resnumber(p_atom))")
#     end
#     mutated_hbond_res = Set{String}()
#     for (p_atom, _, _, _) in mutated_result[:hbonds]
#         push!(mutated_hbond_res, "$(resname(p_atom)) $(resnumber(p_atom))")
#     end
#     native_hbond_str = join(native_hbond_res, ", ")
#     mutated_hbond_str = join(mutated_hbond_res, ", ")

#     native_pi_pi_res = Set{String}()
#     for (p_atom, _, _, _) in native_result[:pi_pi]
#         push!(native_pi_pi_res, "$(resname(p_atom)) $(resnumber(p_atom))")
#     end
#     mutated_pi_pi_res = Set{String}()
#     for (p_atom, _, _, _) in mutated_result[:pi_pi]
#         push!(mutated_pi_pi_res, "$(resname(p_atom)) $(resnumber(p_atom))")
#     end
#     native_pi_pi_str = join(native_pi_pi_res, ", ")
#     mutated_pi_pi_str = join(mutated_pi_pi_res, ", ")

#     native_hbond_pct = native_total > 0 ? round((native_hbonds / native_total) * 100, digits=1) : 0.0
#     mutated_hbond_pct = mutated_total > 0 ? round((mutated_hbonds / mutated_total) * 100, digits=1) : 0.0
#     native_nonbonded_pct = native_total > 0 ? round((native_nonbonded / native_total) * 100, digits=1) : 0.0
#     mutated_nonbonded_pct = mutated_total > 0 ? round((mutated_nonbonded / mutated_total) * 100, digits=1) : 0.0
#     native_insight = native_hbond_pct > native_nonbonded_pct ? "H-bonds dominate ($native_hbond_pct%)" : "Hydrophobic dominate ($native_nonbonded_pct%)"
#     mutated_insight = mutated_hbond_pct > mutated_nonbonded_pct ? "H-bonds dominate ($mutated_hbond_pct%)" : "Hydrophobic dominate ($mutated_nonbonded_pct%)"

#     df = DataFrame(
#         Metric = [
#             "Ligand",
#             "Total Interactions",
#             "Interaction Difference",
#             "Hydrogen Bonds",
#             "H-Bond Difference",
#             "Non-Bonded (Hydrophobic)",
#             "Non-Bonded Difference",
#             "π-π Stacking",
#             "π-π Difference",
#             "Top Interacting Residues",
#             "Changed Residues",
#             "H-Bond Residues",
#             "π-π Residues",
#             "Binding Insight"
#         ],
#         Native = [
#             native_ligand,
#             native_total,
#             "",
#             native_hbonds,
#             "",
#             native_nonbonded,
#             "",
#             native_pi_pi,
#             "",
#             native_top_str,
#             "",
#             native_hbond_str,
#             native_pi_pi_str,
#             native_insight
#         ],
#         Mutated = [
#             mutated_ligand,
#             mutated_total,
#             total_diff,
#             mutated_hbonds,
#             hbond_diff,
#             mutated_nonbonded,
#             nonbonded_diff,
#             mutated_pi_pi,
#             pi_pi_diff,
#             mutated_top_str,
#             changed_residues_str,
#             mutated_hbond_str,
#             mutated_pi_pi_str,
#             mutated_insight
#         ]
#     )

#     CSV.write(output_file, df, missingstring="")
#     println("Side-by-side comparison saved to $output_file")

#     detail_rows = []
#     for (complex_name, result) in complex_results
#         for (p_atom, l_atom, dist, bond_type) in result[:close_contacts]
#             push!(detail_rows, (
#                 Complex = complex_name,
#                 Interaction_Type = "Close Contact",
#                 Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
#                 Protein_Atom = atomname(p_atom),
#                 Ligand_Residue = resname(l_atom),
#                 Ligand_Atom = atomname(l_atom),
#                 Distance = round(dist, digits=2),
#                 Bond_Description = bond_type
#             ))
#         end
#         for (p_atom, l_atom, dist, bond_type) in result[:hbonds]
#             push!(detail_rows, (
#                 Complex = complex_name,
#                 Interaction_Type = "Hydrogen Bond",
#                 Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
#                 Protein_Atom = atomname(p_atom),
#                 Ligand_Residue = resname(l_atom),
#                 Ligand_Atom = atomname(l_atom),
#                 Distance = round(dist, digits=2),
#                 Bond_Description = bond_type
#             ))
#         end
#         for (p_atom, l_atom, dist, bond_type) in result[:nonbonded]
#             push!(detail_rows, (
#                 Complex = complex_name,
#                 Interaction_Type = "Non-Bonded",
#                 Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
#                 Protein_Atom = atomname(p_atom),
#                 Ligand_Residue = resname(l_atom),
#                 Ligand_Atom = atomname(l_atom),
#                 Distance = round(dist, digits=2),
#                 Bond_Description = bond_type
#             ))
#         end
#         for (p_atom, l_atom, dist, bond_type) in result[:pi_pi]
#             push!(detail_rows, (
#                 Complex = complex_name,
#                 Interaction_Type = "π-π Stacking",
#                 Protein_Residue = "$(resname(p_atom)) $(resnumber(p_atom))",
#                 Protein_Atom = atomname(p_atom),
#                 Ligand_Residue = resname(l_atom),
#                 Ligand_Atom = atomname(l_atom),
#                 Distance = round(dist, digits=2),
#                 Bond_Description = bond_type
#             ))
#         end
#     end
#     detail_df = DataFrame(detail_rows)
#     CSV.write(joinpath("public", "outputs", "detailed_interactions.csv"), detail_df)
#     println("Detailed interactions saved to public/outputs/detailed_interactions.csv")

#     residues = sort(collect(all_residues))
#     native_counts = [get(native_result[:interacting_residues], res, 0) for res in residues]
#     mutated_counts = [get(mutated_result[:interacting_residues], res, 0) for res in residues]

#     n = length(residues)
#     x = 1:n
#     bar_width = 0.35
#     plt = bar(
#         x .- bar_width/2, native_counts,
#         label="First Complex ($native_name)",
#         color=:blue,
#         alpha=0.6,
#         bar_width=bar_width,
#         xticks=(1:n, residues),
#         xrotation=45,
#         xlabel="Residue",
#         ylabel="Interaction Count",
#         title="Residue Interaction Counts",
#         legend=:topright,
#         size=(800, 600)
#     )
#     bar!(plt,
#         x .+ bar_width/2, mutated_counts,
#         label="Second Complex ($mutated_name)",
#         color=:red,
#         alpha=0.6,
#         bar_width=bar_width
#     )

#     savefig(plt, plot_file)
#     println("Bar chart saved to $plot_file")
# end

# function print_analytical_summary(complex_results)
#     for (complex_name, result) in complex_results
#         println("\nAnalytical Summary for $complex_name (Ligand: $(result[:ligand_resname])):")
#         println("Total Interactions: ", length(result[:close_contacts]))
#         println("Hydrogen Bonds: ", length(result[:hbonds]))
#         println("Non-Bonded (Hydrophobic): ", length(result[:nonbonded]))
#         println("π-π Stacking: ", length(result[:pi_pi]))
#         sorted_residues = sort(collect(result[:interacting_residues]), by=x->x[2], rev=true)
#         top_residues = sorted_residues[1:min(3, length(sorted_residues))]
#         println("Top Interacting Residues: ", join(["$res ($count)" for (res, count) in top_residues], ", "))
#         hbond_residues = Set{String}()
#         for (p_atom, _, _, _) in result[:hbonds]
#             push!(hbond_residues, "$(resname(p_atom)) $(resnumber(p_atom))")
#         end
#         println("Residues with H-Bonds: ", join(hbond_residues, ", "))
#         pi_pi_residues = Set{String}()
#         for (p_atom, _, _, _) in result[:pi_pi]
#             push!(pi_pi_residues, "$(resname(p_atom)) $(resnumber(p_atom))")
#         end
#         println("Residues with π-π Stacking: ", join(pi_pi_residues, ", "))
#         hbond_percent = length(result[:close_contacts]) > 0 ? round((length(result[:hbonds]) / length(result[:close_contacts])) * 100, digits=1) : 0.0
#         nonbonded_percent = length(result[:close_contacts]) > 0 ? round((length(result[:nonbonded]) / length(result[:close_contacts])) * 100, digits=1) : 0.0
#         println("Binding Insight: ", hbond_percent > nonbonded_percent ? "Hydrogen bonds dominate ($hbond_percent%)" : "Hydrophobic interactions dominate ($nonbonded_percent%)")
#     end
# end

function print_analytical_summary(complex_results)
    ordered_complexes = ["first_complex.cif", "second_complex.cif"]
    for complex_name in ordered_complexes
        if haskey(complex_results, complex_name)
            result = complex_results[complex_name]
            display_name = complex_name == "first_complex.cif" ? "First Complex" : "Second Complex"
            println("\nAnalytical Summary for $display_name (Ligand: $(result[:ligand_resname])):")
            println("Total Interactions: ", length(result[:close_contacts]))
            println("Hydrogen Bonds: ", length(result[:hbonds]))
            println("Non-Bonded (Hydrophobic): ", length(result[:nonbonded]))
            println("π-π Stacking: ", length(result[:pi_pi]))
            sorted_residues = sort(collect(result[:interacting_residues]), by=x->x[2], rev=true)
            top_residues = sorted_residues[1:min(3, length(sorted_residues))]
            println("Top Interacting Residues: ", join(["$res ($count)" for (res, count) in top_residues], ", "))
            hbond_residues = Set{String}()
            for (p_atom, _, _, _) in result[:hbonds]
                push!(hbond_residues, "$(resname(p_atom)) $(resnumber(p_atom))")
            end
            println("Residues with H-Bonds: ", join(hbond_residues, ", "))
            pi_pi_residues = Set{String}()
            for (p_atom, _, _, _) in result[:pi_pi]
                push!(pi_pi_residues, "$(resname(p_atom)) $(resnumber(p_atom))")
            end
            println("Residues with π-π Stacking: ", join(pi_pi_residues, ", "))
            hbond_percent = length(result[:close_contacts]) > 0 ? round((length(result[:hbonds]) / length(result[:close_contacts])) * 100, digits=1) : 0.0
            nonbonded_percent = length(result[:close_contacts]) > 0 ? round((length(result[:nonbonded]) / length(result[:close_contacts])) * 100, digits=1) : 0.0
            println("Binding Insight: ", hbond_percent > nonbonded_percent ? "Hydrogen bonds dominate ($hbond_percent%)" : "Hydrophobic interactions dominate ($nonbonded_percent%)")
        end
    end
end