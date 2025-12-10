# lib/output_results.jl
using CSV
using DataFrames
using CairoMakie
using Printf
using ColorSchemes
using Colors

ENV["GKSwstype"] = "100"

function compare_and_save_results(complex_results, first_key, second_key, comparison_csv_path, detailed_csv_path, plot_path)
    
    all_details = DataFrame()
    for (key, data) in complex_results
        if haskey(data, :hbonds) && !isempty(data[:hbonds])
            df = DataFrame(data[:hbonds])
            df[!, :Complex_Source] .= key
            append!(all_details, df)
        end
    end

    if isempty(all_details)
        CSV.write(detailed_csv_path, DataFrame(Complex_Source=[], Type=[], Residue=[], Distance=[]))
        f = Figure(size = (800, 400)); ax = Axis(f[1, 1]); hidedecorations!(ax); hidespines!(ax)
        text!(ax, 0.5, 0.5, text = "No Interactions Detected", align = (:center, :center), fontsize=16)
        save(plot_path, f)
        return
    end

    desired_cols = [:Complex_Source, :type, :prot_resn, :prot_resi, :prot_atom_name, :lig_atom_name, :distance]
    select!(all_details, intersect(desired_cols, propertynames(all_details)))
    CSV.write(detailed_csv_path, all_details)

    function format_residue(resn, resi); return string(uppercase(string(resn)), resi); end
    for (key, data) in complex_results
        if haskey(data, :interacting_residues)
            new_dict = Dict{String, Int}(); for (res_key, count) in data[:interacting_residues]; parts = split(res_key); if length(parts) == 2; new_dict[format_residue(parts[1], parts[2])] = count; end; end; data[:interacting_residues] = new_dict
        end
    end

    res_counts_1 = complex_results[first_key][:interacting_residues]; res_counts_2 = complex_results[second_key][:interacting_residues]
    all_residues = union(keys(res_counts_1), keys(res_counts_2))
    
    comp_df = DataFrame(Residue = String[], C1 = Int[], C2 = Int[], Total = Int[])
    for res in all_residues; c1 = get(res_counts_1, res, 0); c2 = get(res_counts_2, res, 0); push!(comp_df, (res, c1, c2, c1 + c2)); end
    sort!(comp_df, :Total, rev=true); CSV.write(comparison_csv_path, comp_df) 
    
    plot_df = nrow(comp_df) > 15 ? first(comp_df, 15) : comp_df
    plot_df = sort(plot_df, :Total) 
    
    target_residues = plot_df.Residue
    all_types = sort(unique(all_details[!, :type]))

    calc_height = max(450, 100 + (length(target_residues) * 50))
    f = Figure(size = (900, calc_height), backgroundcolor = :white)
    
    ax = Axis(f[2, 1], 
        title = "Interactome Comparison", titlesize = 18,
        xlabel = "Number of Interactions", xlabelsize = 14,
        yticks = (1:length(target_residues), target_residues), yticklabelsize = 12,
        ygridvisible = false, xgridstyle = :dash, xgridcolor = :gray85
    )

    hidespines!(ax) 

    colors = [colorant"#007bff", colorant"#ff8552", colorant"#28a745", colorant"#dc3545", colorant"#6f42c1", colorant"#17a2b8", colorant"#6c757d"]
    type_color_map = Dict(t => colors[mod1(i, length(colors))] for (i, t) in enumerate(all_types))

    for (i, res) in enumerate(target_residues)
        subset1 = filter(row -> row.Complex_Source == first_key && format_residue(row.prot_resn, row.prot_resi) == res, all_details)
        curr_bot = 0.0
        for type in all_types; cnt = count(row -> row.type == type, eachrow(subset1)); if cnt > 0; barplot!(ax, [i - 0.2], [cnt], offset=[curr_bot], direction=:x, color=type_color_map[type], width=0.35, gap=0.0); curr_bot += cnt; end; end
        subset2 = filter(row -> row.Complex_Source == second_key && format_residue(row.prot_resn, row.prot_resi) == res, all_details)
        curr_bot = 0.0
        for type in all_types; cnt = count(row -> row.type == type, eachrow(subset2)); if cnt > 0; barplot!(ax, [i + 0.2], [cnt], offset=[curr_bot], direction=:x, color=type_color_map[type], width=0.35, gap=0.0); curr_bot += cnt; end; end
    end

    max_x = isempty(plot_df.Total) ? 5 : maximum(plot_df.Total)
    ax.xticks = 0:Int(ceil(max_x))
    CairoMakie.xlims!(ax, 0, max_x + 0.5)
    
    # MOVED THIS BLOCK HERE
    elems = [PolyElement(color = type_color_map[t]) for t in all_types]
    Legend(f[1, 1], elems, all_types, "Interaction Type", nbanks = 4, framevisible = false, labelsize=12, tellheight=true, padding=(0,0,5,0))
    
    Label(f[3, 1], "Upper Bar: Second Complex | Lower Bar: First Complex", fontsize = 10, color = :gray50, justification = :center)
    rowgap!(f.layout, 5)
    resize_to_layout!(f)

    save(plot_path, f)
end

function print_analytical_summary(io::IO, complex_results, first_key, second_key)
    # (Same summary logic)
end