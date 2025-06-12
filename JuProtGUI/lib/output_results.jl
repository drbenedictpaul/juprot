using CSV
using DataFrames
using Plots
using BioStructures

function check_residue_differences(first_structure, second_structure)
    first_model = first_structure[1]
    second_model = second_structure[1]
    first_residues = Dict{String, String}()
    second_residues = Dict{String, String}()

    for res in collectresidues(first_model, standardselector)
        key = "$(resnumber(res))"
        first_residues[key] = resname(res)
    end
    for res in collectresidues(second_model, standardselector)
        key = "$(resnumber(res))"
        second_residues[key] = resname(res)
    end

    differences = String[]
    for key in union(keys(first_residues), keys(second_residues))
        first_res = get(first_residues, key, "None")
        second_res = get(second_residues, key, "None")
        if first_res != second_res
            push!(differences, "Position $key: First=$first_res, Second=$second_res")
        end
    end
    return differences
end

function compare_and_save_results(complex_results, output_file, plot_file, first_structure, second_structure)
    if length(complex_results) != 2
        println("Error: Exactly two complexes required for comparison.")
        return
    end

    residue_diffs = check_residue_differences(first_structure, second_structure)
    if isempty(residue_diffs)
        println("No residue differences detected between complexes.")
    else
        println("Residue differences between complexes:")
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

    total_diff = second_total - first_total
    hbond_diff = second_hbonds - first_hbonds
    nonbonded_diff = second_nonbonded - first_nonbonded

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

    total_specific = first_hbonds + first_nonbonded
    first_hbond_pct = total_specific > 0 ? round((first_hbonds / total_specific) * 100, digits=1) : 0.0
    first_nonbonded_pct = total_specific > 0 ? round((first_nonbonded / total_specific) * 100, digits=1) : 0.0
    total_specific = second_hbonds + second_nonbonded
    second_hbond_pct = total_specific > 0 ? round((second_hbonds / total_specific) * 100, digits=1) : 0.0
    second_nonbonded_pct = total_specific > 0 ? round((second_nonbonded / total_specific) * 100, digits=1) : 0.0
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
            "Top Interacting Residues",
            "Changed Residues",
            "H-Bond Residues",
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
            first_top_str,
            "",
            first_hbond_str,
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
            second_top_str,
            changed_residues_str,
            second_hbond_str,
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
    label="Native or First Complex",
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
    label="Mutated or Second Complex",
    color=:red,
    alpha=0.6,
    bar_width=bar_width
)

    savefig(plt, plot_file)
    println("Bar chart saved to $plot_file")
end

function print_analytical_summary(complex_results)
    ordered_complexes = ["first_complex.cif", "second_complex.cif"]
    for complex_name in ordered_complexes
        if haskey(complex_results, complex_name)
            result = complex_results[complex_name]
            display_name = complex_name == "first_complex.cif" ? "Native or First Complex" : "Mutated or Second Complex"
            total_interactions = length(result[:close_contacts])
            println("\nAnalytical Summary for $display_name (Ligand: $(result[:ligand_resname])):")
            println("Total Interactions: ", total_interactions, " (all atom pairs within 4.0 Å)")
            println("Hydrogen Bonds: ", length(result[:hbonds]))
            println("Non-Bonded (Hydrophobic): ", length(result[:nonbonded]))
            sorted_residues = sort(collect(result[:interacting_residues]), by=x->x[2], rev=true)
            top_residues = sorted_residues[1:min(3, length(sorted_residues))]
            println("Top Interacting Residues: ", join(["$res ($count)" for (res, count) in top_residues], ", "))
            hbond_residues = Set{String}()
            for (p_atom, _, _, _) in result[:hbonds]
                push!(hbond_residues, "$(resname(p_atom)) $(resnumber(p_atom))")
            end
            println("Residues with H-Bonds: ", join(hbond_residues, ", "))
            hbond_percent = total_interactions > 0 ? round((length(result[:hbonds]) / total_interactions) * 100, digits=1) : 0.0
            nonbonded_percent = total_interactions > 0 ? round((length(result[:nonbonded]) / total_interactions) * 100, digits=1) : 0.0
            println("Binding Insight: ", hbond_percent > nonbonded_percent ? "Hydrogen bonds dominate ($hbond_percent%)" : "Hydrophobic interactions dominate ($nonbonded_percent%)")
        end
    end
end