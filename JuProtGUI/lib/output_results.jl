# lib/output_results.jl
using CSV
using DataFrames
using Plots

"""
    compare_and_save_results(complex_results, first_complex_key, second_complex_key, output_csv_path, detailed_csv_path, plot_file_path)

Compares H-bond interaction results from two complexes and saves them to CSV and a plot.
"""
function compare_and_save_results(
        complex_results::Dict{String, Dict}, 
        first_complex_key::String, 
        second_complex_key::String, 
        output_csv_path::String, 
        detailed_csv_path::String,
        plot_file_path::String
    )

    if !haskey(complex_results, first_complex_key) || !haskey(complex_results, second_complex_key)
        err_msg = "Error: One or both complex keys ('$first_complex_key', '$second_complex_key') not found in results."
        println(err_msg)
        # Create empty/error files to prevent downstream serving errors
        open(io -> println(io, err_msg), output_csv_path, "w")
        open(io -> println(io, err_msg), detailed_csv_path, "w")
        plt_err = plot(title="Error: Input Data Missing", size=(800,600))
        annotate!(plt_err, 0.5, 0.5, Plots.text("Could not generate plot due to missing complex data.", 10, :center))
        savefig(plt_err, plot_file_path)
        return # Exit the function
    end

    # Check if results themselves contain error messages (from process_pdb_files loop)
    if haskey(complex_results[first_complex_key], "error_message") || haskey(complex_results[second_complex_key], "error_message")
        err_msg1 = get(complex_results[first_complex_key], "error_message", "")
        err_msg2 = get(complex_results[second_complex_key], "error_message", "")
        full_err_msg =strip("Error processing complexes: $err_msg1; $err_msg2", [';',' '])
        println(full_err_msg)
        open(io -> println(io, full_err_msg), output_csv_path, "w")
        open(io -> println(io, full_err_msg), detailed_csv_path, "w")
        plt_err = plot(title="Error: PDB Processing Failed", size=(800,600))
        annotate!(plt_err, 0.5, 0.5, Plots.text(full_err_msg, 8, :center))
        savefig(plt_err, plot_file_path)
        return
    end


    first_result = complex_results[first_complex_key]
    second_result = complex_results[second_complex_key]

    first_ligand = first_result[:ligand_resname]
    second_ligand = second_result[:ligand_resname]

    first_total_hbonds = length(first_result[:hbonds])
    second_total_hbonds = length(second_result[:hbonds])
    hbond_diff = second_total_hbonds - first_total_hbonds

    sort_dict_by_value(d) = isempty(d) ? [] : sort(collect(d), by=x->x[2], rev=true)
    
    first_sorted_interactions = sort_dict_by_value(first_result[:interacting_residues])
    second_sorted_interactions = sort_dict_by_value(second_result[:interacting_residues])
    
    first_top_str = isempty(first_sorted_interactions) ? "None" : join(["$(res) ($(count))" for (res, count) in first_sorted_interactions[1:min(3, end)]], ", ")
    second_top_str = isempty(second_sorted_interactions) ? "None" : join(["$(res) ($(count))" for (res, count) in second_sorted_interactions[1:min(3, end)]], ", ")

    all_interacting_residues_keys = union(keys(first_result[:interacting_residues]), keys(second_result[:interacting_residues]))
    changed_residues_info = String[]
    for res_key in all_interacting_residues_keys
        first_count = get(first_result[:interacting_residues], res_key, 0)
        second_count = get(second_result[:interacting_residues], res_key, 0)
        if first_count != second_count
            push!(changed_residues_info, "$res_key: First=$first_count, Second=$second_count")
        end
    end
    changed_residues_str = isempty(changed_residues_info) ? "None" : join(changed_residues_info, "; ")

    get_hbond_residue_set(bonds) = Set{String}(["$(bond.prot_resn) $(bond.prot_resi)" for bond in bonds])
    
    first_hbond_res_set = get_hbond_residue_set(first_result[:hbonds])
    second_hbond_res_set = get_hbond_residue_set(second_result[:hbonds])
    
    first_hbond_str = isempty(first_hbond_res_set) ? "None" : join(sort(collect(first_hbond_res_set)), ", ")
    second_hbond_str = isempty(second_hbond_res_set) ? "None" : join(sort(collect(second_hbond_res_set)), ", ")

    df_summary = DataFrame(
        Metric = [
            "Ligand Name",
            "Total H-Bonds Detected",
            "Difference in H-Bonds (Second - First)",
            "Top Interacting Residues (ResID Count)",
            "Residues with Changed Interaction Counts",
            "Unique Residues Forming H-Bonds"
        ],
        First_Complex = [
            first_ligand,
            first_total_hbonds,
            "", 
            first_top_str,
            "", 
            first_hbond_str
        ],
        Second_Complex = [
            second_ligand,
            second_total_hbonds,
            hbond_diff,
            second_top_str,
            changed_residues_str,
            second_hbond_str
        ]
    )

    CSV.write(output_csv_path, df_summary, missingstring="N/A")
    println("Side-by-side comparison summary saved to $output_csv_path")

    detail_rows = []
    complex_data_map = [
        (name="First Complex ($first_ligand)", data=first_result),
        (name="Second Complex ($second_ligand)", data=second_result)
    ]

    for item in complex_data_map
        complex_display_name = item.name
        result_data = item.data
        
        for bond_detail in result_data[:hbonds]
            push!(detail_rows, (
                Complex_Name = complex_display_name,
                Interaction_Type = bond_detail.type,
                Protein_Residue_ID = "$(bond_detail.prot_resn) $(bond_detail.prot_resi)",
                Protein_Atom_PDB_Name = bond_detail.prot_atom_name,
                Protein_Atom_Type = bond_detail.prot_atom_type,
                Ligand_Residue_Name = result_data[:ligand_resname],
                Ligand_Atom_PDB_Name = bond_detail.lig_atom_name,
                Ligand_Atom_Type = bond_detail.lig_atom_type,
                Distance_Angstrom = round(bond_detail.distance, digits=2),
                Angle_Degrees = round(bond_detail.angle, digits=2)
            ))
        end
    end
    
    if !isempty(detail_rows)
        detail_df = DataFrame(detail_rows)
        CSV.write(detailed_csv_path, detail_df)
        println("Detailed H-bond interactions saved to $detailed_csv_path")
    else
        empty_df_headers = DataFrame(
            Complex_Name=String[], Interaction_Type=String[], Protein_Residue_ID=String[],
            Protein_Atom_PDB_Name=String[], Protein_Atom_Type=String[], Ligand_Residue_Name=String[],
            Ligand_Atom_PDB_Name=String[], Ligand_Atom_Type=String[], Distance_Angstrom=Float64[], Angle_Degrees=Float64[]
        )
        CSV.write(detailed_csv_path, empty_df_headers)
        println("No detailed H-bond interactions to save; empty CSV created at $detailed_csv_path")
    end

    sorted_residue_keys = sort(collect(all_interacting_residues_keys))

    if !isempty(sorted_residue_keys)
        first_counts = [get(first_result[:interacting_residues], res_key, 0) for res_key in sorted_residue_keys]
        second_counts = [get(second_result[:interacting_residues], res_key, 0) for res_key in sorted_residue_keys]

        num_residues = length(sorted_residue_keys)
        x_ticks_positions = 1:num_residues
        
        bar_width = num_residues > 20 ? 0.4 : 0.35
        fig_width = max(800, num_residues * 40) 
        fig_height = 600

        plt = bar(
            x_ticks_positions .- bar_width/2, first_counts,
            label="First Complex ($first_ligand)",
            color=:dodgerblue, 
            alpha=0.7,
            bar_width=bar_width,
            xticks=(x_ticks_positions, sorted_residue_keys),
            xrotation=60, 
            xlabel="Protein Residue (ResName ResID)",
            ylabel="Number of H-Bonds",
            title="Comparison of H-Bond Counts per Residue",
            legend=:topright,
            size=(fig_width, fig_height),
            left_margin=50Plots.px, 
            bottom_margin=100Plots.px 
        )
        bar!(plt,
            x_ticks_positions .+ bar_width/2, second_counts,
            label="Second Complex ($second_ligand)",
            color=:orangered, 
            alpha=0.7,
            bar_width=bar_width
        )
        all_counts = vcat(first_counts, second_counts)
        ylims!(0, max(1, isempty(all_counts) ? 1 : maximum(all_counts)) * 1.1)

        savefig(plt, plot_file_path)
        println("Bar chart of residue H-bond counts saved to $plot_file_path")
    else
        plt = plot(title="No Interacting Residues Found", size=(800,600))
        annotate!(plt, 0.5, 0.5, Plots.text("No H-bond interactions were found for either complex.", 10, :center))
        savefig(plt, plot_file_path)
        println("No interacting residues to plot; placeholder image saved to $plot_file_path")
    end
end

"""
    print_analytical_summary(io::IO, complex_results, first_complex_key, second_complex_key)

Prints an analytical summary of H-bond interactions to the given IO stream.
"""
function print_analytical_summary(io::IO, complex_results::Dict{String, Dict}, first_complex_key::String, second_complex_key::String)
    println(io, "Analytical Summary of Protein-Ligand Hydrogen Bonds")
    println(io, "===================================================")

    complex_keys_to_print = [first_complex_key, second_complex_key]
    
    for (idx, complex_key) in enumerate(complex_keys_to_print)
        if haskey(complex_results, complex_key)
            result = complex_results[complex_key]
            # Check if this result itself is an error marker from process_pdb_files
            if haskey(result, "error_message")
                display_name_err = idx == 1 ? "First Complex ($first_complex_key)" : "Second Complex ($second_complex_key)"
                println(io, "\n--- $display_name_err ---")
                println(io, "Error during processing: ", result["error_message"])
                continue # Skip to next complex
            end

            display_name = idx == 1 ? "First Complex" : "Second Complex"
            
            println(io, "\n--- $display_name (Ligand: $(result[:ligand_resname])) ---")
            println(io, "Original PDB source: $(basename(get(result, :original_path, "N/A")))")
            
            num_hbonds = length(result[:hbonds])
            println(io, "Total Hydrogen Bonds Detected: ", num_hbonds)

            if num_hbonds > 0
                sorted_interacting_residues = isempty(result[:interacting_residues]) ? [] : sort(collect(result[:interacting_residues]), by=x->x[2], rev=true)
                
                top_residues_str_parts = String[]
                for res_info in sorted_interacting_residues[1:min(5, end)]
                     push!(top_residues_str_parts, "$(res_info[1]) (involved in $(res_info[2]) H-bond(s))")
                end
                top_residues_str = join(top_residues_str_parts, "\n  - ")


                println(io, "Top Interacting Residues (by H-bond count):")
                if !isempty(top_residues_str) print(io, "  - "); println(io, top_residues_str) else println(io, "  None") end


                hbond_forming_residues = Set{String}()
                for bond_detail in result[:hbonds]
                    push!(hbond_forming_residues, "$(bond_detail.prot_resn) $(bond_detail.prot_resi)")
                end
                sorted_hbond_residues = sort(collect(hbond_forming_residues))
                println(io, "Unique Residues Forming H-Bonds: ", isempty(sorted_hbond_residues) ? "None" : join(sorted_hbond_residues, ", "))
                
                println(io, "Sample H-Bond Details (up to 3):")
                if isempty(result[:hbonds])
                    println(io, "  None")
                else
                    for bond in result[:hbonds][1:min(3, end)]
                        println(io, "  - Protein: $(bond.prot_atom_name) ($(bond.prot_resn) $(bond.prot_resi)) <--> Ligand: $(bond.lig_atom_name) ($(result[:ligand_resname])) | Dist: $(round(bond.distance,digits=2))Å, Angle: $(round(bond.angle,digits=2))°")
                    end
                end
            else
                println(io, "No hydrogen bonds were identified for this complex with the specified ligand.")
            end
        else
            display_name_err = idx == 1 ? "First Complex ($first_complex_key)" : "Second Complex ($second_complex_key)"
            println(io, "\nWarning: Data for $display_name_err not found in results.")
        end
    end
    println(io, "\n===================================================")
    println(io, "Summary generated by juProt.")
end