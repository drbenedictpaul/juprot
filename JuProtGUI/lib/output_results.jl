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
        println("Error: Exactly two complexes (native and mutated) required for comparison.")
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

    complexes = collect(keys(complex_results))
    native_name, mutated_name = complexes[1], complexes[2]
    native_result = complex_results[native_name]
    mutated_result = complex_results[mutated_name]

    native_ligand = native_result[:ligand_resname]
    mutated_ligand = mutated_result[:ligand_resname]

    native_total = length(native_result[:close_contacts])
    mutated_total = length(mutated_result[:close_contacts])
    native_hbonds = length(native_result[:hbonds])
    mutated_hbonds = length(mutated_result[:hbonds])
    native_nonbonded = length(native_result[:nonbonded])
    mutated_nonbonded = length(mutated_result[:nonbonded])
    native_pi_pi = length(native_result[:pi_pi])
    mutated_pi_pi = length(mutated_result[:pi_pi])

    total_diff = mutated_total - native_total
    hbond_diff = mutated_hbonds - native_hbonds
    nonbonded_diff = mutated_nonbonded - native_nonbonded
    pi_pi_diff = mutated_pi_pi - native_pi_pi

    native_sorted = sort(collect(native_result[:interacting_residues]), by=x->x[2], rev=true)
    mutated_sorted = sort(collect(mutated_result[:interacting_residues]), by=x->x[2], rev=true)
    native_top = native_sorted[1:min(3, length(native_sorted))]
    mutated_top = mutated_sorted[1:min(3, length(mutated_sorted))]
    native_top_str = join(["$res ($count)" for (res, count) in native_top], ", ")
    mutated_top_str = join(["$res ($count)" for (res, count) in mutated_top], ", ")

    all_residues = union(keys(native_result[:interacting_residues]), keys(mutated_result[:interacting_residues]))
    changed_residues = String[]
    for res in all_residues
        native_count = get(native_result[:interacting_residues], res, 0)
        mutated_count = get(mutated_result[:interacting_residues], res, 0)
        if native_count != mutated_count
            push!(changed_residues, "$res: Native=$native_count, Mutated=$mutated_count")
        end
    end
    changed_residues_str = join(changed_residues, "; ")

    native_hbond_res = Set{String}()
    for (p_atom, _, _, _) in native_result[:hbonds]
        push!(native_hbond_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    mutated_hbond_res = Set{String}()
    for (p_atom, _, _, _) in mutated_result[:hbonds]
        push!(mutated_hbond_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    native_hbond_str = join(native_hbond_res, ", ")
    mutated_hbond_str = join(mutated_hbond_res, ", ")

    native_pi_pi_res = Set{String}()
    for (p_atom, _, _, _) in native_result[:pi_pi]
        push!(native_pi_pi_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    mutated_pi_pi_res = Set{String}()
    for (p_atom, _, _, _) in mutated_result[:pi_pi]
        push!(mutated_pi_pi_res, "$(resname(p_atom)) $(resnumber(p_atom))")
    end
    native_pi_pi_str = join(native_pi_pi_res, ", ")
    mutated_pi_pi_str = join(mutated_pi_pi_res, ", ")

    native_hbond_pct = native_total > 0 ? round((native_hbonds / native_total) * 100, digits=1) : 0.0
    mutated_hbond_pct = mutated_total > 0 ? round((mutated_hbonds / mutated_total) * 100, digits=1) : 0.0
    native_nonbonded_pct = native_total > 0 ? round((native_nonbonded / native_total) * 100, digits=1) : 0.0
    mutated_nonbonded_pct = mutated_total > 0 ? round((mutated_nonbonded / mutated_total) * 100, digits=1) : 0.0
    native_insight = native_hbond_pct > native_nonbonded_pct ? "H-bonds dominate ($native_hbond_pct%)" : "Hydrophobic dominate ($native_nonbonded_pct%)"
    mutated_insight = mutated_hbond_pct > mutated_nonbonded_pct ? "H-bonds dominate ($mutated_hbond_pct%)" : "Hydrophobic dominate ($mutated_nonbonded_pct%)"

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
        Native = [
            native_ligand,
            native_total,
            "",
            native_hbonds,
            "",
            native_nonbonded,
            "",
            native_pi_pi,
            "",
            native_top_str,
            "",
            native_hbond_str,
            native_pi_pi_str,
            native_insight
        ],
        Mutated = [
            mutated_ligand,
            mutated_total,
            total_diff,
            mutated_hbonds,
            hbond_diff,
            mutated_nonbonded,
            nonbonded_diff,
            mutated_pi_pi,
            pi_pi_diff,
            mutated_top_str,
            changed_residues_str,
            mutated_hbond_str,
            mutated_pi_pi_str,
            mutated_insight
        ]
    )

    CSV.write(output_file, df, missingstring="")
    println("Side-by-side comparison saved to $output_file")

    detail_rows = []
    for (complex_name, result) in complex_results
        for (p_atom, l_atom, dist, bond_type) in result[:close_contacts]
            push!(detail_rows, (
                Complex = complex_name,
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
                Complex = complex_name,
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
                Complex = complex_name,
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
                Complex = complex_name,
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
    CSV.write("detailed_interactions.csv", detail_df)
    println("Detailed interactions saved to detailed_interactions.csv")

    residues = sort(collect(all_residues))
    native_counts = [get(native_result[:interacting_residues], res, 0) for res in residues]
    mutated_counts = [get(mutated_result[:interacting_residues], res, 0) for res in residues]

    n = length(residues)
    x = 1:n
    bar_width = 0.35
    plt = bar(
        x .- bar_width/2, native_counts,
        label="Native ($native_name)",
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
        x .+ bar_width/2, mutated_counts,
        label="Mutated ($mutated_name)",
        color=:red,
        alpha=0.6,
        bar_width=bar_width
    )

    savefig(plt, plot_file)
    println("Bar chart saved to $plot_file")
end

function print_analytical_summary(complex_results)
    for (complex_name, result) in complex_results
        println("\nAnalytical Summary for $complex_name (Ligand: $(result[:ligand_resname])):")
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