# using CSV
# using DataFrames
# using Plots
# using Printf

# # --- 1. COMPARISON & SAVING ---
# function compare_and_save_results(complex_results, first_key, second_key, comparison_csv_path, detailed_csv_path, plot_path)
    
#     all_details = DataFrame()
#     for (key, data) in complex_results
#         if haskey(data, :hbonds) && !isempty(data[:hbonds])
#             df = DataFrame(data[:hbonds])
#             df[!, :Complex_Source] .= key
#             append!(all_details, df)
#         end
#     end

#     if isempty(all_details)
#         CSV.write(detailed_csv_path, DataFrame(Complex_Source=[],Type=[],Residue=[]))
#         return
#     end

#     desired_cols = [:Complex_Source, :type, :prot_resn, :prot_resi, :prot_chain]
#     select!(all_details, intersect(desired_cols, propertynames(all_details)))
#     CSV.write(detailed_csv_path, all_details)

#     # Comparison Table
#     res1 = complex_results[first_key][:interacting_residues]
#     res2 = complex_results[second_key][:interacting_residues]
#     all_res = union(keys(res1), keys(res2))
    
#     cdf = DataFrame(Residue=String[], C1=Int[], C2=Int[], Total=Int[])
#     for res in all_res
#         c1 = get(res1, res, 0)
#         c2 = get(res2, res, 0)
#         push!(cdf, (res, c1, c2, c1+c2))
#     end
#     sort!(cdf, :Total, rev=true)
#     CSV.write(comparison_csv_path, cdf)
# end

# # --- 2. STACKED CHART ---
# function generate_stacked_chart(cr, fk, sk, output_path)
#     all_residues = Set{String}()
#     function get_res_data(key)
#         counts = Dict{String, Dict{String, Int}}()
#         for d in cr[key][:hbonds]
#             c_display = isempty(d.prot_chain) ? "_" : d.prot_chain
#             rk = "$(c_display) $(d.prot_resn) $(d.prot_resi)"
#             push!(all_residues, rk)
#             t = titlecase(replace(string(d.type), "_" => " "))
#             if !haskey(counts, rk); counts[rk] = Dict{String, Int}(); end
#             counts[rk][t] = get(counts[rk], t, 0) + 1
#         end
#         return counts
#     end

#     data1 = get_res_data(fk); data2 = get_res_data(sk)
    
#     # Sort residues naturally
#     sorted_res = sort(collect(all_residues), by=x -> tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
    
#     if isempty(sorted_res)
#         p = plot(title="No Interactions", framestyle=:box); savefig(p, output_path); return
#     end

#     int_types = ["Hydrogen Bond", "Hydrophobic", "Salt Bridge", "Pi Stacking", "Pi Cation", "Halogen", "Water Bridge", "Metal"]
#     type_colors = Dict("Hydrogen Bond" => :blue, "Hydrophobic" => :orange, "Salt Bridge" => :red, "Pi Stacking" => :green, "Pi Cation" => :lightgreen, "Halogen" => :purple, "Water Bridge"=> :cyan, "Metal" => :gray, "Other" => :black)

#     function build_matrix(data_dict)
#         y_vals = Dict{String, Vector{Int}}()
#         for t in int_types
#             y_vals[t] = [get(get(data_dict, r, Dict()), t, 0) for r in sorted_res]
#         end
#         return y_vals
#     end

#     y1 = build_matrix(data1); y2 = build_matrix(data2)

#     function make_subplot(y_data, title_str)
#         p = plot(title=title_str, xrotation=45, legend=:outertopright, ylabel="Count")
#         active_types = [t for t in int_types if sum(y_data[t]) > 0]
#         if isempty(active_types)
#             bar!(p, sorted_res, zeros(length(sorted_res)), label="")
#         else
#             mat = hcat([y_data[t] for t in active_types]...)
#             cols = reshape([type_colors[t] for t in active_types], 1, :)
#             labels = reshape(active_types, 1, :)
#             bar!(p, sorted_res, mat, label=labels, color=cols)
#         end
#         return p
#     end

#     p1 = make_subplot(y1, "Complex 1: $(cr[fk][:ligand_resname])")
#     p2 = make_subplot(y2, "Complex 2: $(cr[sk][:ligand_resname])")
#     final_plot = plot(p1, p2, layout=(2,1), size=(1000, 800), margin=10Plots.mm)
#     savefig(final_plot, output_path)
# end

# # --- 3. ANALYTICAL SUMMARY (Robust Version) ---
# function print_analytical_summary(io, cr, fk, sk)
#     r1 = cr[fk]; r2 = cr[sk]
#     l1 = get(r1, :hbonds, []); l2 = get(r2, :hbonds, [])
#     c1 = length(l1); c2 = length(l2)
    
#     println(io, "=== Comparative Analysis Summary ===")
#     println(io, "Total Interactions: $c1 (First) vs $c2 (Second)")
#     if c1 > 0 && c2 > 0
#         println(io, "Fold Change: $(round(c2/c1, digits=2))-fold (Complex 2 vs Complex 1)")
#     end

#     # Breakdown by Type
#     function count_types(list)
#         d = Dict{String,Int}()
#         for i in list; d[i.type] = get(d, i.type, 0) + 1; end
#         return d
#     end
#     t1 = count_types(l1); t2 = count_types(l2)
#     println(io, "\n--- Breakdown by Interaction Type ---")
#     println(io, rpad("Type", 20) * " | " * rpad("Complex 1", 10) * " | " * rpad("Complex 2", 10))
#     println(io, "-"^48)
#     for t in sort(collect(union(keys(t1), keys(t2))))
#         println(io, rpad(t, 20) * " | " * rpad(string(get(t1, t, 0)), 10) * " | " * rpad(string(get(t2, t, 0)), 10))
#     end

#     # Residue Differences (Robust Normalization)
#     println(io, "\n--- Significant Residue Differences ---")
    
#     function normalize_residues(res_dict)
#         norm = Dict{String, Int}()
#         for (k, count) in res_dict
#             # k is usually "CHAIN RESN RESI" (e.g. "_ ARG 55" or "A ARG 55")
#             # split() removes empty strings, so "_ ARG 55" -> ["_", "ARG", "55"]
#             parts = split(k)
#             if length(parts) >= 2
#                 # Always take the last two parts: RESN and RESI
#                 # This ignores the chain ID at the start
#                 clean_key = "$(parts[end-1]) $(parts[end])"
#                 norm[clean_key] = get(norm, clean_key, 0) + count
#             end
#         end
#         return norm
#     end

#     n_res1 = normalize_residues(r1[:interacting_residues])
#     n_res2 = normalize_residues(r2[:interacting_residues])
    
#     all_keys = sort(collect(union(keys(n_res1), keys(n_res2))), by=x->tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
    
#     changes = false
#     for res in all_keys
#         v1 = get(n_res1, res, 0); v2 = get(n_res2, res, 0)
#         if v1 != v2
#             changes = true
#             if v1 == 0; println(io, "  [GAINED]  $res: 0 -> $v2")
#             elseif v2 == 0; println(io, "  [LOST]    $res: $v1 -> 0")
#             else; println(io, "  [CHANGED] $res: $v1 -> $v2"); end
#         end
#     end
#     if !changes; println(io, "  No residue-level changes observed."); end
# end

# # --- 4. 2D SPATIAL MAP (Fuzzy Search Version) ---
# function generate_2d_network_map(cr, fk, sk, output_path)
    
#     function get_geometry(pdb_path, lig_name, lig_meta, interaction_list)
#         target_residues = Dict{Int, String}()
#         for i in interaction_list
#             target_residues[i.prot_resi] = isempty(i.prot_chain) ? "_" : i.prot_chain
#         end
        
#         lig_target_resi = get(lig_meta, "resi", 0)
#         lig_coords = []; res_coords = Dict{String, Tuple{Float64, Float64}}()
        
#         for line in eachline(pdb_path)
#             if startswith(line, "ENDMDL"); break; end # Stop at end of first model
#             if startswith(line, "ATOM") || startswith(line, "HETATM")
#                 try
#                     # Loose parsing to handle various formats
#                     line_len = length(line)
#                     if line_len < 46; continue; end # Skip short lines
                    
#                     resi_str = line[23:26]
#                     resi = tryparse(Int, resi_str)
                    
#                     atom_name = strip(line[13:16])
#                     res_name  = strip(line[17:20])
                    
#                     # Chain extraction
#                     chain = "_"
#                     if line_len >= 22 && !isspace(line[22]); chain = string(line[22]); end
                    
#                     x = parse(Float64, line[31:38])
#                     y = parse(Float64, line[39:46])
                    
#                     # --- LIGAND CAPTURE (Fuzzy) ---
#                     # 1. Match by Residue Number (Strongest)
#                     is_ligand = false
#                     if resi == lig_target_resi
#                          is_ligand = true
#                     # 2. Fallback: Match by Name if ResNumber failed (e.g. UNK vs LIG)
#                     elseif occursin(lig_name, line) || (lig_name != "UNK" && occursin(lig_name, res_name))
#                          is_ligand = true
#                     end
                    
#                     if is_ligand
#                         push!(lig_coords, (x, y))
#                     end
                    
#                     # --- PROTEIN RESIDUE CAPTURE ---
#                     if resi !== nothing && haskey(target_residues, resi)
#                         req_chain = target_residues[resi]
#                         chain_match = (req_chain == "_") || (chain == "_") || (req_chain == chain)
                        
#                         if chain_match && atom_name == "CA"
#                             label = "$(res_name)$(resi)"
#                             res_coords[label] = (x, y)
#                         end
#                     end
#                 catch; continue; end
#             end
#         end
#         return (lig_coords, res_coords)
#     end

#     function create_subplot(key, title_text)
#         data = cr[key]
#         (lig_points, res_points) = get_geometry(data[:original_path], data[:ligand_resname], data[:ligand_data], data[:hbonds])
        
#         plt = plot(title=title_text, framestyle=:box, grid=true, aspect_ratio=:equal, xlabel="X (Å)", ylabel="Y (Å)", legend=:topright)
        
#         lx, ly = 0.0, 0.0
#         if !isempty(lig_points)
#             lx = sum(p[1] for p in lig_points) / length(lig_points)
#             ly = sum(p[2] for p in lig_points) / length(lig_points)
#             scatter!(plt, [lx], [ly], marker=:star5, markersize=12, color=:yellow, markerstrokecolor=:black, label="Ligand")
#         else
#              # Fallback: Draw center at 0,0 if ligand not found
#              scatter!(plt, [0], [0], marker=:star5, markersize=12, color=:gray, label="Ligand (Center N/A)")
#              annotate!(plt, 0, 1, Plots.text("Ligand Geometry Missing", 8, :red))
#         end
            
#         if !isempty(res_points)
#             for (label, coords) in res_points
#                 rx, ry = coords
#                 res_num = parse(Int, filter(isdigit, label))
#                 my_int = [x for x in data[:hbonds] if x.prot_resi == res_num]
#                 itype = isempty(my_int) ? "Other" : my_int[1].type
                
#                 col = occursin("Hydrogen", itype) ? :blue : 
#                       occursin("Hydrophobic", itype) ? :orange :
#                       occursin("Salt", itype) ? :red :
#                       occursin("Pi", itype) ? :green : :gray
                
#                 plot!(plt, [lx, rx], [ly, ry], line=(1, :dash, col), label="")
#                 scatter!(plt, [rx], [ry], marker=:circle, markersize=8, color=col, label="")
                
#                 # Label Offset
#                 dx, dy = rx - lx, ry - ly
#                 len = sqrt(dx^2 + dy^2)
#                 if len > 0; dx, dy = dx/len, dy/len; end
#                 annotate!(plt, rx + dx*2, ry + dy*2, Plots.text(label, 8, :black))
#             end
#         end
#         return plt
#     end

#     p1 = create_subplot(fk, "Complex 1")
#     p2 = create_subplot(sk, "Complex 2")
#     l = plot(framestyle=:none, grid=false, ticks=false, showaxis=false)
#     plot!(l, [0], [0], label="H-Bond", color=:blue, lw=2)
#     plot!(l, [0], [0], label="Hydrophobic", color=:orange, lw=2)
#     plot!(l, [0], [0], label="Salt Bridge", color=:red, lw=2)
#     plot!(l, [0], [0], label="Pi-Stacking", color=:green, lw=2)
    
#     final = plot(p1, p2, l, layout=@layout([a b; c{0.1h}]), size=(900, 600), margin=5Plots.mm)
#     savefig(final, output_path)
# end

# # --- 5. HTML TABLES ---
# function generate_interaction_tables(cr, fk, sk)
#     function make_table(key)
#         data = cr[key]
#         res_map = Dict{String, Vector{String}}()
#         for d in data[:hbonds]
#             c_display = isempty(d.prot_chain) ? "_" : d.prot_chain
#             rk = "$(c_display) $(d.prot_resn) $(d.prot_resi)"
#             t = titlecase(replace(string(d.type), "_" => " "))
#             if !haskey(res_map, rk); res_map[rk] = String[]; end
#             push!(res_map[rk], t)
#         end
#         rows = ""
#         sorted_keys = sort(collect(keys(res_map)), by=x->tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
#         for rk in sorted_keys
#             types = res_map[rk]
#             counts = Dict{String, Int}()
#             for t in types; counts[t] = get(counts, t, 0) + 1; end
#             desc = join(["$(v) $(k)" for (k,v) in counts], ", ")
#             rows *= "<tr><td style='border:1px solid #ddd;padding:6px;'>$rk</td><td style='border:1px solid #ddd;padding:6px;'>$desc</td></tr>"
#         end
#         return "<table style='width:100%;border-collapse:collapse;font-size:0.9em;'><thead><tr style='background:#f8f9fa;'><th style='border:1px solid #ddd;padding:6px;text-align:left;'>Residue</th><th style='border:1px solid #ddd;padding:6px;text-align:left;'>Interaction Types</th></tr></thead><tbody>$rows</tbody></table>"
#     end
#     t1 = make_table(fk); t2 = make_table(sk)
#     return """<div style="display:flex; gap:20px; flex-wrap:wrap;"><div style="flex:1; min-width:300px;"><h3 style="color:#007bff; border-bottom:2px solid #eee; padding-bottom:10px;">First Complex ($(cr[fk][:ligand_resname]))</h3>$t1</div><div style="flex:1; min-width:300px;"><h3 style="color:#007bff; border-bottom:2px solid #eee; padding-bottom:10px;">Second Complex ($(cr[sk][:ligand_resname]))</h3>$t2</div></div>"""
# end

# Filename: lib/output_results.jl

using CSV
using DataFrames
using Plots
using Printf

# --- 1. COMPARISON & SAVING ---
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
        CSV.write(detailed_csv_path, DataFrame(Complex_Source=[],Type=[],Residue=[]))
        return
    end

    desired_cols = [:Complex_Source, :type, :prot_resn, :prot_resi, :prot_chain]
    select!(all_details, intersect(desired_cols, propertynames(all_details)))
    CSV.write(detailed_csv_path, all_details)

    res1 = complex_results[first_key][:interacting_residues]
    res2 = complex_results[second_key][:interacting_residues]
    all_res = union(keys(res1), keys(res2))
    
    cdf = DataFrame(Residue=String[], C1=Int[], C2=Int[], Total=Int[])
    for res in all_res
        c1 = get(res1, res, 0)
        c2 = get(res2, res, 0)
        push!(cdf, (res, c1, c2, c1+c2))
    end
    sort!(cdf, :Total, rev=true)
    CSV.write(comparison_csv_path, cdf)
end

# --- 2. STACKED CHART ---
function generate_stacked_chart(cr, fk, sk, output_path)
    all_residues = Set{String}()
    function get_res_data(key)
        counts = Dict{String, Dict{String, Int}}()
        for d in cr[key][:hbonds]
            c_display = isempty(d.prot_chain) ? "_" : d.prot_chain
            rk = "$(c_display) $(d.prot_resn) $(d.prot_resi)"
            push!(all_residues, rk)
            t = titlecase(replace(string(d.type), "_" => " "))
            if !haskey(counts, rk); counts[rk] = Dict{String, Int}(); end
            counts[rk][t] = get(counts[rk], t, 0) + 1
        end
        return counts
    end

    data1 = get_res_data(fk); data2 = get_res_data(sk)
    sorted_res = sort(collect(all_residues), by=x -> tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
    
    if isempty(sorted_res)
        p = plot(title="No Interactions", framestyle=:box); savefig(p, output_path); return
    end

    int_types = ["Hydrogen Bond", "Hydrophobic", "Salt Bridge", "Pi Stacking", "Pi Cation", "Halogen", "Water Bridge", "Metal"]
    type_colors = Dict("Hydrogen Bond" => :blue, "Hydrophobic" => :orange, "Salt Bridge" => :red, "Pi Stacking" => :green, "Pi Cation" => :lightgreen, "Halogen" => :purple, "Water Bridge"=> :cyan, "Metal" => :gray, "Other" => :black)

    function build_matrix(data_dict)
        y_vals = Dict{String, Vector{Int}}()
        for t in int_types; y_vals[t] = [get(get(data_dict, r, Dict()), t, 0) for r in sorted_res]; end
        return y_vals
    end

    y1 = build_matrix(data1); y2 = build_matrix(data2)

    function make_subplot(y_data, title_str)
        p = plot(title=title_str, xrotation=45, legend=:outertopright, ylabel="Count")
        active_types = [t for t in int_types if sum(y_data[t]) > 0]
        if isempty(active_types)
            bar!(p, sorted_res, zeros(length(sorted_res)), label="")
        else
            mat = hcat([y_data[t] for t in active_types]...)
            cols = reshape([type_colors[t] for t in active_types], 1, :)
            labels = reshape(active_types, 1, :)
            bar!(p, sorted_res, mat, label=labels, color=cols)
        end
        return p
    end

    p1 = make_subplot(y1, "Complex 1: $(cr[fk][:ligand_resname])")
    p2 = make_subplot(y2, "Complex 2: $(cr[sk][:ligand_resname])")
    final_plot = plot(p1, p2, layout=(2,1), size=(1000, 800), margin=10Plots.mm)
    savefig(final_plot, output_path)
end

# --- 3. ANALYTICAL SUMMARY ---
function print_analytical_summary(io, cr, fk, sk)
    r1 = cr[fk]; r2 = cr[sk]
    l1 = get(r1, :hbonds, []); l2 = get(r2, :hbonds, [])
    c1 = length(l1); c2 = length(l2)
    
    println(io, "=== Comparative Analysis Summary ===")
    println(io, "Total Interactions: $c1 (First) vs $c2 (Second)")
    if c1 > 0 && c2 > 0; println(io, "Fold Change: $(round(c2/c1, digits=2))-fold (C2 vs C1)"); end

    count_types(list) = begin; d=Dict{String,Int}(); for i in list; d[i.type]=get(d,i.type,0)+1; end; return d; end
    t1 = count_types(l1); t2 = count_types(l2)
    println(io, "\n--- Breakdown by Interaction Type ---")
    println(io, rpad("Type", 20) * " | " * rpad("Complex 1", 10) * " | " * "Complex 2")
    println(io, "-"^48)
    for t in sort(collect(union(keys(t1), keys(t2)))); println(io, rpad(t, 20) * " | " * rpad(string(get(t1,t,0)), 10) * " | " * string(get(t2,t,0))); end

    println(io, "\n--- Significant Residue Differences ---")
    normalize_residues(rd) = begin; n=Dict{String,Int}(); for(k,c) in rd; p=split(k); if length(p)>=2; nk="$(p[end-1]) $(p[end])"; n[nk]=get(n,nk,0)+c; end; end; return n; end
    n1=normalize_residues(r1[:interacting_residues]); n2=normalize_residues(r2[:interacting_residues])
    all_keys = sort(collect(union(keys(n1), keys(n2))), by=x->tryparse(Int,split(x)[end])===nothing ? 0 : parse(Int,split(x)[end]))
    
    changes=false
    for res in all_keys
        v1=get(n1,res,0); v2=get(n2,res,0)
        if v1!=v2; changes=true; if v1==0;println(io,"  [GAINED]  $res: 0 -> $v2");elseif v2==0;println(io,"  [LOST]    $res: $v1 -> 0");else;println(io,"  [CHANGED] $res: $v1 -> $v2");end;end
    end
    if !changes; println(io, "  No residue-level changes observed."); end
end

# --- 4. 2D SPATIAL MAP ---
function generate_2d_network_map(cr, fk, sk, output_path)
    function get_geometry(pdb_path, lig_name, lig_meta, interaction_list)
        target_residues = Dict{Int, String}()
        for i in interaction_list; target_residues[i.prot_resi] = isempty(i.prot_chain) ? "_" : i.prot_chain; end
        lig_target_resi = get(lig_meta, "resi", 0)
        lig_coords=[]; res_coords=Dict{String,Tuple{Float64,Float64}}()
        
        for line in eachline(pdb_path)
            if startswith(line,"ENDMDL"); break; end
            if startswith(line,"ATOM")||startswith(line,"HETATM")
                try
                    resi=tryparse(Int,line[23:26]); atom_name=strip(line[13:16]); res_name=strip(line[17:20])
                    chain = length(line)>=22 && !isspace(line[22]) ? string(line[22]) : "_"
                    x=parse(Float64,line[31:38]); y=parse(Float64,line[39:46])
                    
                    if (resi==lig_target_resi) || (occursin(lig_name,line)); push!(lig_coords,(x,y)); end
                    
                    if resi!==nothing && haskey(target_residues,resi) && (target_residues[resi]=="_" || chain=="_" || target_residues[resi]==chain) && atom_name=="CA"
                        res_coords["$(res_name)$(resi)"]=(x,y)
                    end
                catch; continue; end
            end
        end
        return (lig_coords,res_coords)
    end

    function create_subplot(key, title_text)
        data = cr[key]
        (lig_points, res_points) = get_geometry(data[:original_path], data[:ligand_resname], data[:ligand_data], data[:hbonds])
        plt = plot(title=title_text, framestyle=:box, aspect_ratio=:equal, xlabel="X (Å)", ylabel="Y (Å)", legend=:topright)
        
        lx,ly = isempty(lig_points) ? (0.0,0.0) : (sum(p[1] for p in lig_points)/length(lig_points), sum(p[2] for p in lig_points)/length(lig_points))
        scatter!(plt, [lx], [ly], marker=:star5, markersize=12, color=isempty(lig_points) ? :gray : :yellow, markerstrokecolor=:black, label="Ligand")
        
        for (label, coords) in res_points
            rx,ry = coords; res_num = parse(Int,filter(isdigit,label)); my_int = [x for x in data[:hbonds] if x.prot_resi==res_num]
            itype = isempty(my_int) ? "Other" : my_int[1].type
            col = occursin("Hydrogen",itype) ? :blue : occursin("Hydrophobic",itype) ? :orange : occursin("Salt",itype) ? :red : occursin("Pi",itype) ? :green : :gray
            plot!(plt,[lx,rx],[ly,ry],line=(1,:dash,col),label=""); scatter!(plt,[rx],[ry],marker=:circle,markersize=8,color=col,label="")
            dx,dy=rx-lx,ry-ly; len=sqrt(dx^2+dy^2); if len>0;dx,dy=dx/len,dy/len;end; annotate!(plt,rx+dx*2,ry+dy*2,Plots.text(label,8,:black))
        end
        return plt
    end

    p1 = create_subplot(fk, "Complex 1"); p2 = create_subplot(sk, "Complex 2")
    l = plot(framestyle=:none,grid=false,ticks=false,showaxis=false)
    plot!(l,[0],[0],label="H-Bond",color=:blue,lw=2); plot!(l,[0],[0],label="Hydrophobic",color=:orange,lw=2)
    plot!(l,[0],[0],label="Salt Bridge",color=:red,lw=2); plot!(l,[0],[0],label="Pi-Stacking",color=:green,lw=2)
    
    final = plot(p1, p2, l, layout=@layout([a b; c{0.1h}]), size=(900, 600), margin=5Plots.mm)
    savefig(final, output_path)
end

# --- 5. HTML TABLES ---
function generate_interaction_tables(cr, fk, sk)
    function make_table(key)
        res_map = Dict{String,Vector{String}}()
        for d in cr[key][:hbonds]
            rk = "$(isempty(d.prot_chain) ? "_" : d.prot_chain) $(d.prot_resn) $(d.prot_resi)"
            t = titlecase(replace(string(d.type), "_" => " "))
            if !haskey(res_map, rk); res_map[rk] = []; end; push!(res_map[rk], t)
        end
        rows = ""
        sorted_keys = sort(collect(keys(res_map)), by=x->tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
        for rk in sorted_keys
            counts=Dict{String,Int}(); for t in res_map[rk]; counts[t]=get(counts,t,0)+1; end
            desc = join(["$(v) $(k)" for(k,v) in counts], ", ")
            rows *= "<tr><td style='border:1px solid #ddd;padding:6px;'>$rk</td><td style='border:1px solid #ddd;padding:6px;'>$desc</td></tr>"
        end
        return "<table style='width:100%;border-collapse:collapse;font-size:0.9em;'><thead><tr style='background:#f8f9fa;'><th style='border:1px solid #ddd;padding:6px;text-align:left;'>Residue</th><th style='border:1px solid #ddd;padding:6px;text-align:left;'>Interactions</th></tr></thead><tbody>$rows</tbody></table>"
    end
    t1 = make_table(fk); t2 = make_table(sk)
    return """<div style="display:flex;gap:20px;flex-wrap:wrap;"><div style="flex:1;min-width:300px;"><h3 style="color:#007bff;border-bottom:2px solid #eee;padding-bottom:10px;">Complex 1 ($(cr[fk][:ligand_resname]))</h3>$t1</div><div style="flex:1;min-width:300px;"><h3 style="color:#007bff;border-bottom:2px solid #eee;padding-bottom:10px;">Complex 2 ($(cr[sk][:ligand_resname]))</h3>$t2</div></div>"""
end