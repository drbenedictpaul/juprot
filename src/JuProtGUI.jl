# module JuProtGUI
#     # 1. FORCE SILENT MODE (Prevents Cloud Run Crash)
#     ENV["GKSwstype"] = "100"

#     using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
#     using CSV, DataFrames, PythonCall, Printf, Downloads, Plots, Base64

#     # --- CONFIGURATION ---
#     const LIB_PATH = joinpath(@__DIR__, "..", "lib")
    
#     # We use /tmp for temporary generation
#     const OUTPUT_DIR = "/tmp/outputs"
#     mkpath(OUTPUT_DIR)
    
#     include(joinpath(LIB_PATH, "output_results.jl"))
#     include(joinpath(LIB_PATH, "detect_hbonds.jl"))

#     # --- STACKED CHART GENERATOR ---
#     function generate_stacked_chart(cr, fk, sk, output_path)
#         all_residues = Set{String}()
#         function get_res_data(key)
#             counts = Dict{String, Dict{String, Int}}()
#             for d in cr[key][:hbonds]
#                 rk = "$(d.prot_resn) $(d.prot_resi)"
#                 push!(all_residues, rk)
#                 t_raw = hasproperty(d, :type) ? string(d.type) : "Other"
#                 t = titlecase(replace(t_raw, "_" => " "))
#                 if !haskey(counts, rk); counts[rk] = Dict{String, Int}(); end
#                 counts[rk][t] = get(counts[rk], t, 0) + 1
#             end
#             return counts
#         end

#         data1 = get_res_data(fk); data2 = get_res_data(sk)
#         sorted_res = sort(collect(all_residues), by=x->tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
        
#         if isempty(sorted_res)
#             p = plot(title="No Interactions Found", framestyle=:box)
#             savefig(p, output_path)
#             return
#         end

#         int_types = ["Hydrogen Bond", "Hydrophobic", "Salt Bridge", "Pi Stacking", "Pi Cation", "Halogen", "Water Bridge", "Metal"]
#         type_colors = Dict("Hydrogen Bond" => :blue, "Hydrophobic" => :orange, "Salt Bridge" => :red, "Pi Stacking" => :green, "Pi Cation" => :lightgreen, "Halogen" => :purple, "Water Bridge"=> :cyan, "Metal" => :gray, "Other" => :black)

#         function build_matrix(data_dict)
#             y_vals = Dict{String, Vector{Int}}()
#             for t in int_types
#                 y_vals[t] = [get(get(data_dict, r, Dict()), t, 0) for r in sorted_res]
#             end
#             return y_vals
#         end

#         y1 = build_matrix(data1); y2 = build_matrix(data2)

#         function make_subplot(y_data, title_str)
#             p = plot(title=title_str, xrotation=45, legend=:outertopright, ylabel="Count")
#             active_types = [t for t in int_types if sum(y_data[t]) > 0]
#             if isempty(active_types)
#                 bar!(p, sorted_res, zeros(length(sorted_res)), label="")
#             else
#                 mat = hcat([y_data[t] for t in active_types]...)
#                 cols = reshape([type_colors[t] for t in active_types], 1, :)
#                 labels = reshape(active_types, 1, :)
#                 # bar!(p, sorted_res, mat, label=labels, color=cols, bar_position=:stack)
#                 bar!(p, sorted_res, mat, label=labels, color=cols)
#             end
#             return p
#         end

#         p1 = make_subplot(y1, "Complex 1: $(cr[fk][:ligand_resname])")
#         p2 = make_subplot(y2, "Complex 2: $(cr[sk][:ligand_resname])")
#         final_plot = plot(p1, p2, layout=(2,1), size=(1000, 800), margin=5Plots.mm)
#         savefig(final_plot, output_path)
#     end

#     # --- HTML TABLE GENERATOR ---
#     function generate_interaction_tables(cr, fk, sk)
#         function make_table(key)
#             data = cr[key]
#             res_map = Dict{String, Vector{String}}()
#             for d in data[:hbonds]
#                 rk = "$(d.prot_resn) $(d.prot_resi)"
#                 t_raw = hasproperty(d, :type) ? string(d.type) : "Interaction"
#                 t = titlecase(replace(t_raw, "_" => " "))
#                 if !haskey(res_map, rk); res_map[rk] = String[]; end
#                 push!(res_map[rk], t)
#             end
#             rows = ""
#             sorted_keys = sort(collect(keys(res_map)), by=x->tryparse(Int, split(x)[end]) === nothing ? 0 : parse(Int, split(x)[end]))
#             for rk in sorted_keys
#                 types = res_map[rk]
#                 counts = Dict{String, Int}()
#                 for t in types; counts[t] = get(counts, t, 0) + 1; end
#                 desc = join(["$(v) $(k)" for (k,v) in counts], ", ")
#                 rows *= "<tr><td style='border:1px solid #ddd;padding:6px;'>$rk</td><td style='border:1px solid #ddd;padding:6px;'>$desc</td></tr>"
#             end
#             return "<table style='width:100%;border-collapse:collapse;font-size:0.9em;'><thead><tr style='background:#f8f9fa;'><th style='border:1px solid #ddd;padding:6px;text-align:left;'>Residue</th><th style='border:1px solid #ddd;padding:6px;text-align:left;'>Interaction Types</th></tr></thead><tbody>$rows</tbody></table>"
#         end
#         t1 = make_table(fk); t2 = make_table(sk)
#         return """<div style="display:flex; gap:20px; flex-wrap:wrap;"><div style="flex:1; min-width:300px;"><h3 style="color:#007bff; border-bottom:2px solid #eee; padding-bottom:10px;">First Complex ($(cr[fk][:ligand_resname]))</h3>$t1</div><div style="flex:1; min-width:300px;"><h3 style="color:#007bff; border-bottom:2px solid #eee; padding-bottom:10px;">Second Complex ($(cr[sk][:ligand_resname]))</h3>$t2</div></div>"""
#     end

#     function generate_pymol_script(complex_results)
#         script_content = IOBuffer()
        
#         path1 = complex_results["first_complex.pdb"][:original_path]
#         path2 = complex_results["second_complex.pdb"][:original_path]
        
#         # 1. Retrieve Ligand Data
#         l1_data = complex_results["first_complex.pdb"][:ligand_data]
#         l2_data = complex_results["second_complex.pdb"][:ligand_data]
        
#         # 2. Retrieve Ligand Names (Needed for the fix)
#         ln1 = complex_results["first_complex.pdb"][:ligand_resname]
#         ln2 = complex_results["second_complex.pdb"][:ligand_resname]
        
#         name1 = "Complex1"
#         name2 = "Complex2"

#         println(script_content, "reinitialize")
#         println(script_content, "bg_color white")
        
#         # --- Python Block for Safe Loading ---
#         println(script_content, "python")
#         println(script_content, "import pymol")
#         println(script_content, "from pymol import cmd")
        
#         if isfile(path1)
#             println(script_content, "pdb1_string = \"\"\"")
#             for line in eachline(path1); println(script_content, line); end
#             println(script_content, "\"\"\"")
#             println(script_content, "cmd.read_pdbstr(pdb1_string, '$name1')")
#         end

#         if isfile(path2)
#             println(script_content, "pdb2_string = \"\"\"")
#             for line in eachline(path2); println(script_content, line); end
#             println(script_content, "\"\"\"")
#             println(script_content, "cmd.read_pdbstr(pdb2_string, '$name2')")
#         end
        
#         println(script_content, "python end")

#         # --- Visualization ---
#         println(script_content, "hide everything")
#         println(script_content, "show cartoon")
#         println(script_content, "color gray80")
#         println(script_content, "align $name2, $name1")

#         # Helper for Protein Residues (Keep Chain logic here as proteins usually have chains)
#         function get_selection_string(d)
#             sel_parts = String[]
#             for k in keys(d)
#                 parts = split(k)
#                 if length(parts) >= 3
#                     c_id = parts[1]
#                     r_id = parts[3]
#                     if c_id == "_"
#                         push!(sel_parts, "(resi $r_id)")
#                     else
#                         push!(sel_parts, "(chain $c_id and resi $r_id)")
#                     end
#                 end
#             end
#             if isempty(sel_parts); return "none"; end
#             return join(sel_parts, " or ")
#         end

#         r1_sel = get_selection_string(complex_results["first_complex.pdb"][:interacting_residues])
#         r2_sel = get_selection_string(complex_results["second_complex.pdb"][:interacting_residues])

#         if r1_sel != "none"
#             println(script_content, "select int_res1, ($name1 and ($r1_sel))")
#             println(script_content, "show sticks, int_res1")
#             println(script_content, "color marine, int_res1")
#             println(script_content, "util.cnc int_res1")
#         end

#         if r2_sel != "none"
#             println(script_content, "select int_res2, ($name2 and ($r2_sel))")
#             println(script_content, "show sticks, int_res2")
#             println(script_content, "color firebrick, int_res2")
#             println(script_content, "util.cnc int_res2")
#         end

#         # --- FIXED LIGAND SELECTION ---
#         # We rely on Residue Number AND Residue Name.
#         # We deliberately IGNORE the Chain ID to avoid mismatches in docked files.
#         function make_lig_sel(name, data, lig_name)
#             r = data["resi"]
#             return "($name and resi $r and resn $lig_name)"
#         end

#         l1_sel_cmd = make_lig_sel(name1, l1_data, ln1)
#         l2_sel_cmd = make_lig_sel(name2, l2_data, ln2)

#         println(script_content, "select lig1, $l1_sel_cmd")
#         println(script_content, "select lig2, $l2_sel_cmd")
        
#         println(script_content, "show sticks, lig1")
#         println(script_content, "color marine, lig1")
#         println(script_content, "util.cnc lig1")

#         println(script_content, "show sticks, lig2")
#         println(script_content, "color firebrick, lig2")
#         println(script_content, "util.cnc lig2")

#         println(script_content, "zoom lig1 or lig2, 5")
#         println(script_content, "deselect")
        
#         return String(take!(script_content))
#     end
    
#     function get_ligand_names(p)
#         try
#             if !isfile(p); return ["ERROR: File missing at $p"]; end
#             plip=pyimport("plip.structure.preparation"); mol=plip.PDBComplex(); mol.load_pdb(p); l=Set{String}()
#             if PythonCall.pyhasattr(mol,"ligands")&&mol.ligands!==PythonCall.pybuiltins.None
#                 for r in mol.ligands
#                     if PythonCall.pyhasattr(r,"hetid"); n_py=r.hetid; if n_py!==PythonCall.pybuiltins.None; n=pyconvert(String,n_py); if !isempty(n)&&n!="HOH"; push!(l,n); end; end; end
#                 end
#             end
#             if isempty(l)&&PythonCall.pyhasattr(mol,"res_het")
#                 for r in mol.res_het; if PythonCall.pyhasattr(r,"hetid"); n=pyconvert(String,r.hetid); if !isempty(n)&&n!="HOH"; push!(l,n); end; end; end
#             end
#             if isempty(l); return ["No Ligands Detected"]; end
#             return collect(l)
#         catch e; err_msg = sprint(showerror, e); return ["ERROR: $err_msg"]; end
#     end
    
#     function process_pdb_files(p1, p2, l1, l2)
#         cr = Dict{String, Dict}()
#         fk = "first_complex.pdb"
#         sk = "second_complex.pdb"
#         files = [(path=p1, ligand=l1, key=fk), (path=p2, ligand=l2, key=sk)]
#         ok = true
        
#         for i in files
#             try
#                 (ints, lig_meta) = detect_hbonds(i.path, i.ligand)
#                 ir = Dict{String, Int}()
#                 for d in ints
#                      c_display = isempty(d.prot_chain) ? "_" : d.prot_chain
#                      rk = "$(c_display) $(d.prot_resn) $(d.prot_resi)"
#                      ir[rk] = get(ir, rk, 0) + 1
#                 end
#                 cr[i.key] = Dict(
#                     :hbonds => ints,
#                     :interacting_residues => ir,
#                     :ligand_resname => i.ligand,
#                     :ligand_data => lig_meta,
#                     :original_path => i.path
#                 )
#             catch e; println("Error processing $(i.path): $e"); ok=false; end
#         end
        
#         if !ok; return Dict("error" => "Processing failed."); end
        
#         session_id = string(abs(rand(Int)))[1:6]
#         ocsv = joinpath(OUTPUT_DIR, "comparison_table_$(session_id).csv")
#         dcsv = joinpath(OUTPUT_DIR, "detailed_interactions_$(session_id).csv")
#         ppath = joinpath(OUTPUT_DIR, "residue_interactions_$(session_id).png")
#         net_path = joinpath(OUTPUT_DIR, "network_map_$(session_id).png") # NEW MAP
 
#         try
#             generate_stacked_chart(cr, fk, sk, ppath)
#             generate_2d_network_map(cr, fk, sk, net_path) # Call the new function
#             compare_and_save_results(cr, fk, sk, ocsv, dcsv, "dummy.png") 
#             st = capture_analytical_summary(cr, fk, sk)
#             tables_html = generate_interaction_tables(cr, fk, sk)
            
#             comp_b64 = base64encode(read(ocsv))
#             det_b64 = base64encode(read(dcsv))
#             plot_b64 = base64encode(read(ppath))
#             net_b64 = base64encode(read(net_path)) # Encode Map
            
#             rm(ocsv, force=true); rm(dcsv, force=true); rm(ppath, force=true); rm(net_path, force=true)
 
#             return Dict(
#                 "comp_b64" => comp_b64, 
#                 "det_b64" => det_b64, 
#                 "plot_b64" => plot_b64, 
#                 "net_b64" => net_b64, 
#                 "summary" => st, 
#                 "tables_html" => tables_html, 
#                 "session_data" => cr
#             )
#         catch e; return Dict("error" => "Visualization failed: $(sprint(showerror, e))"); end
#     end
    
#     function capture_analytical_summary(cr,fk,sk);try;io=IOBuffer();print_analytical_summary(io,cr,fk,sk);return String(take!(io));catch;return"Summary failed.";end;end
 
#     # route("/") do
#     #  html("""
#     #  <!DOCTYPE html><html><head><title>juProt</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.form-group{margin-bottom:20px}label{display:block;margin-bottom:8px;font-weight:700}input[type=file],input[type=text]{width:calc(100% - 22px);padding:10px;border:1px solid #ccc;border-radius:4px}button{padding:10px 20px;background:#007bff;color:#fff;border:none;cursor:pointer;border-radius:4px;font-size:16px}button:hover{background:#0056b3}h1{color:#0056b3;text-align:center}a{color:#007bff;text-decoration:none}.footer-links-container{display:flex;justify-content:center;gap:30px;margin-top:40px;padding-top:20px;border-top:1px solid #eee}.divider{text-align:center;font-weight:700;color:#aaa;margin:20px 0}.validation-popup{display:none;position:absolute;top:35px;left:110px;background-color:#fff;border:1px solid #ccc;border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,.15);padding:8px 12px;z-index:100;white-space:nowrap;font-size:14px;color:#333}.validation-popup::before{content:'';position:absolute;bottom:100%;left:20px;border-width:7px;border-style:solid;border-color:transparent transparent #ccc transparent}.validation-popup::after{content:'';position:absolute;bottom:100%;left:21px;border-width:6px;border-style:solid;border-color:transparent transparent #fff transparent}.validation-popup .icon{display:inline-block;background-color:#ff8552;color:#fff;width:16px;height:16px;border-radius:3px;text-align:center;font-weight:700;line-height:16px;margin-right:8px;font-size:12px}</style><script>function validateForm(){document.getElementById('alert_1').style.display='none';document.getElementById('alert_2').style.display='none';var f1=document.getElementById('first_pdb').value;var i1=document.getElementById('first_pdb_id').value.trim();var c1=(f1!==""||i1!=="");var f2=document.getElementById('second_pdb').value;var i2=document.getElementById('second_pdb_id').value.trim();var c2=(f2!==""||i2!=="");var v=true;if(!c1){document.getElementById('alert_1').style.display='flex';v=false}if(!c2){document.getElementById('alert_2').style.display='flex';v=false}return v}</script></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>juProt: Protein-Ligand Interaction Analyzer</h1><p style="text-align:center;color:#666">Comparative analysis of the complete protein-ligand interactome.</p><form action=/select-ligands method=post enctype=multipart/form-data onsubmit="return validateForm()"><div class=form-group><label for=first_pdb>Upload First Complex (PDB File):</label><input type=file id=first_pdb name=first_pdb accept=.pdb><div id=alert_1 class=validation-popup><span class=icon>!</span>Please select a file or enter a PDB ID.</div></div><div class=form-group><label for=second_pdb>Upload Second Complex (PDB File):</label><input type=file id=second_pdb name=second_pdb accept=.pdb><div id=alert_2 class=validation-popup><span class=icon>!</span>Please select a file or enter a PDB ID.</div></div><div class=divider>OR</div><div class=form-group><label for=first_pdb_id>Enter First PDB ID (e.g., 3EQM):</label><input type=text id=first_pdb_id name=first_pdb_id placeholder="PDB ID"></div><div class=form-group><label for=second_pdb_id>Enter Second PDB ID:</label><input type=text id=second_pdb_id name=second_pdb_id placeholder="PDB ID"></div><button type=submit>Load Ligands</button></form><div class=footer-links-container><p><a href=/how-to-use>How to Use & Applications</a></p><p><a href=/about>About juProt</a></p></div></div></body></html>
#     #  """)
#     # end

#     route("/") do
#      html("""
#      <!DOCTYPE html><html><head><title>juProt</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon">
#      <style>
#         body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}
#         .container{max-width:800px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}
#         .form-group{margin-bottom:20px}
#         label{display:block;margin-bottom:8px;font-weight:700}
#         input[type=file],input[type=text]{width:calc(100% - 22px);padding:10px;border:1px solid #ccc;border-radius:4px}
#         button{padding:10px 20px;background:#007bff;color:#fff;border:none;cursor:pointer;border-radius:4px;font-size:16px}
#         button:hover{background:#0056b3}
#         h1{color:#0056b3;text-align:center}
#         a{color:#007bff;text-decoration:none}
#         .footer-links-container{display:flex;justify-content:center;gap:30px;margin-top:40px;padding-top:20px;border-top:1px solid #eee}
#         .divider{text-align:center;font-weight:700;color:#aaa;margin:20px 0}
#         .validation-popup{display:none;position:absolute;top:35px;left:110px;background-color:#fff;border:1px solid #ccc;border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,.15);padding:8px 12px;z-index:100;white-space:nowrap;font-size:14px;color:#333}.validation-popup::before{content:'';position:absolute;bottom:100%;left:20px;border-width:7px;border-style:solid;border-color:transparent transparent #ccc transparent}.validation-popup::after{content:'';position:absolute;bottom:100%;left:21px;border-width:6px;border-style:solid;border-color:transparent transparent #fff transparent}.validation-popup .icon{display:inline-block;background-color:#ff8552;color:#fff;width:16px;height:16px;border-radius:3px;text-align:center;font-weight:700;line-height:16px;margin-right:8px;font-size:12px}
        
#         /* LOADING OVERLAY CSS */
#         #loader {
#             display: none;
#             position: fixed; top: 0; left: 0; width: 100%; height: 100%;
#             background: rgba(255,255,255,0.9);
#             z-index: 1000;
#             flex-direction: column; justify-content: center; align-items: center;
#         }
#         .spinner {
#             border: 8px solid #f3f3f3; border-top: 8px solid #007bff;
#             border-radius: 50%; width: 60px; height: 60px;
#             animation: spin 1s linear infinite;
#         }
#         @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
#         .loading-text { margin-top: 20px; font-size: 18px; color: #0056b3; font-weight: bold; }
#      </style>
#      <script>
#         function validateForm(){
#             document.getElementById('alert_1').style.display='none';
#             document.getElementById('alert_2').style.display='none';
#             var f1=document.getElementById('first_pdb').value;
#             var i1=document.getElementById('first_pdb_id').value.trim();
#             var c1=(f1!==""||i1!=="");
#             var f2=document.getElementById('second_pdb').value;
#             var i2=document.getElementById('second_pdb_id').value.trim();
#             var c2=(f2!==""||i2!=="");
#             var v=true;
#             if(!c1){document.getElementById('alert_1').style.display='flex';v=false}
#             if(!c2){document.getElementById('alert_2').style.display='flex';v=false}
            
#             // SHOW LOADER IF VALID
#             if(v) {
#                 document.getElementById('loader').style.display = 'flex';
#             }
#             return v;
#         }
#      </script>
#      </head>
#      <body>
#         <!-- LOADER DIV -->
#         <div id="loader">
#             <div class="spinner"></div>
#             <div class="loading-text">Uploading and Parsing Structures...<br><span style="font-size:14px;color:#666;font-weight:normal">This may take a few seconds.</span></div>
#         </div>

#         <div class=container>
#             <div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div>
#             <h1>juProt: Protein-Ligand Interaction Analyzer</h1>
#             <p style="text-align:center;color:#666">Comparative analysis of the complete protein-ligand interactome.</p>
#             <form action=/select-ligands method=post enctype=multipart/form-data onsubmit="return validateForm()">
#                 <div class=form-group>
#                     <label for=first_pdb>Upload First Complex (PDB File):</label>
#                     <input type=file id=first_pdb name=first_pdb accept=.pdb>
#                     <div id=alert_1 class=validation-popup><span class=icon>!</span>Please select a file or enter a PDB ID.</div>
#                 </div>
#                 <div class=form-group>
#                     <label for=second_pdb>Upload Second Complex (PDB File):</label>
#                     <input type=file id=second_pdb name=second_pdb accept=.pdb>
#                     <div id=alert_2 class=validation-popup><span class=icon>!</span>Please select a file or enter a PDB ID.</div>
#                 </div>
#                 <div class=divider>OR</div>
#                 <div class=form-group><label for=first_pdb_id>Enter First PDB ID (e.g., 3EQM):</label><input type=text id=first_pdb_id name=first_pdb_id placeholder="PDB ID"></div>
#                 <div class=form-group><label for=second_pdb_id>Enter Second PDB ID:</label><input type=text id=second_pdb_id name=second_pdb_id placeholder="PDB ID"></div>
#                 <button type=submit>Load Ligands</button>
#             </form>
#             <div class=footer-links-container>
#                 <p><a href=/how-to-use>How to Use & Applications</a></p>
#                 <p><a href=/about>About juProt</a></p>
#             </div>
#         </div>
#      </body></html>
#      """)
#     end
    
#     route("/select-ligands", method=POST) do
#         first_payload=haskey(filespayload(),"first_pdb") ? filespayload("first_pdb") : nothing;
#         second_payload=haskey(filespayload(),"second_pdb") ? filespayload("second_pdb") : nothing
#         f_id=postpayload(:first_pdb_id,""); s_id=postpayload(:second_pdb_id,"")
#         mkpath(OUTPUT_DIR); tp1=""; tp2=""; fn1=""; fn2=""
#         try
#             if !isempty(f_id); tp1=joinpath(OUTPUT_DIR,"$(f_id).pdb"); Downloads.download("https://files.rcsb.org/download/$(f_id).pdb",tp1); fn1="$(f_id).pdb"
#             elseif first_payload!==nothing&&!isempty(first_payload.data); tp1=joinpath(OUTPUT_DIR,"temp_first_"*string(randn())[3:8]*".pdb"); write(tp1,first_payload.data); fn1=first_payload.name
#             else; return html("<h1>Error</h1><p>No input for first complex.</p>"); end
#             if !isempty(s_id); tp2=joinpath(OUTPUT_DIR,"$(s_id).pdb"); Downloads.download("https://files.rcsb.org/download/$(s_id).pdb",tp2); fn2="$(s_id).pdb"
#             elseif second_payload!==nothing&&!isempty(second_payload.data); tp2=joinpath(OUTPUT_DIR,"temp_second_"*string(randn())[3:8]*".pdb"); write(tp2,second_payload.data); fn2=second_payload.name
#             else; return html("<h1>Error</h1><p>No input for second complex.</p>"); end
#             lig1=get_ligand_names(tp1); lig2=get_ligand_names(tp2)
#             if isempty(lig1)&&isempty(lig2); return html("<h1>Error</h1><p>No ligands found.</p>"); end
#             opt1=join(["<option value='$l'>$l</option>" for l in lig1],""); opt2=join(["<option value='$l'>$l</option>" for l in lig2],"")
            
#             # --- UPDATED HTML WITH LOADER ---
#             html_body="""<!DOCTYPE html><html><head><title>juProt: Select Ligands</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon">
#             <style>
#                 body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}
#                 .container{max-width:800px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}
#                 .form-group{margin-bottom:20px}label{display:block;margin-bottom:8px;font-weight:700}
#                 select{width:100%;padding:10px;border:1px solid #ccc;border-radius:4px;box-sizing:border-box}
#                 button{padding:10px 20px;background:#007bff;color:#fff;border:none;cursor:pointer;border-radius:4px;font-size:16px}
#                 button:hover{background:#0056b3}h1{color:#0056b3;text-align:center}.footer-link{text-align:center;margin-top:20px}a{color:#007bff;text-decoration:none}
                
#                 /* LOADER CSS */
#                 #loader {
#                     display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%;
#                     background: rgba(255,255,255,0.9); z-index: 1000;
#                     flex-direction: column; justify-content: center; align-items: center;
#                 }
#                 .spinner {
#                     border: 8px solid #f3f3f3; border-top: 8px solid #007bff;
#                     border-radius: 50%; width: 60px; height: 60px;
#                     animation: spin 1s linear infinite;
#                 }
#                 @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
#                 .loading-text { margin-top: 20px; font-size: 18px; color: #0056b3; font-weight: bold; }
#             </style>
#             <script>
#                 function showLoader() {
#                     document.getElementById('loader').style.display = 'flex';
#                 }
#             </script>
#             </head>
#             <body>
#                 <!-- LOADER DIV -->
#                 <div id="loader">
#                     <div class="spinner"></div>
#                     <div class="loading-text">Calculating Interactions...<br><span style="font-size:14px;color:#666;font-weight:normal">This involves complex calculations (H-bonds, Pi-Stacking, etc).<br>Please wait...</span></div>
#                 </div>

#                 <div class=container>
#                     <div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div>
#                     <h1>Select Target Ligands</h1>
#                     <form action=/analyze method=post onsubmit="showLoader()">
#                         <input type=hidden name=first_pdb_path value="$tp1">
#                         <input type=hidden name=second_pdb_path value="$tp2">
#                         <div class=form-group>
#                             <label>First Complex Ligand ($fn1):</label>
#                             <select name=first_ligand required>$opt1</select>
#                         </div>
#                         <div class=form-group>
#                             <label>Second Complex Ligand ($fn2):</label>
#                             <select name=second_ligand required>$opt2</select>
#                         </div>
#                         <button type=submit>Run Analysis</button>
#                     </form>
#                     <div class=footer-link><p><a href=/>Back</a></p></div>
#                 </div>
#             </body></html>"""
#             return html(html_body)
#         catch e; return html("<h1>Error</h1><p>$(sprint(showerror,e))</p>"); end
#     end
    
#     SESSION_DATA=Ref{Any}(nothing)
#     # route("/analyze",method=POST) do
#     #     p1=postpayload(:first_pdb_path);p2=postpayload(:second_pdb_path);l1=postpayload(:first_ligand);l2=postpayload(:second_ligand)
#     #     result=process_pdb_files(p1,p2,l1,l2)
#     #     if haskey(result,"session_data");SESSION_DATA[]=result["session_data"];end
#     #     if haskey(result,"error");return html("<h1>Error</h1><p>$(result["error"])</p>");end
#     #     ss=replace(result["summary"],"&"=>"&amp;","<"=>"&lt;",">"=>"&gt;","\n"=>"<br>")
#     #     html_body="""<!DOCTYPE html><html><head><title>juProt Analysis Results</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:950px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.result-section{margin-top:20px;padding-bottom:20px;border-bottom:1px solid #eee}.result-section:last-child{border-bottom:none}h1{color:#0056b3;text-align:center;margin-bottom:30px}h2{color:#007bff;margin-bottom:5px;margin-top:0}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}img{max-width:100%;height:auto;border:none;margin-top:0;display:block;margin-left:auto;margin-right:auto}pre{background:#e9ecef;padding:15px;border-radius:5px;white-space:pre-wrap;word-wrap:break-word;font-size:.9em}ul{list-style-type:none;padding-left:0}ul li{margin-bottom:8px}ul li a{display:inline-block;padding:8px 12px;background-color:#007bff;color:#fff;border-radius:4px}ul li a:hover{background-color:#0056b3}.footer-link{text-align:center;margin-top:30px}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>Analysis Results</h1><p>Comparison of full interaction profiles for the selected ligands.</p><div class=result-section><h2>Downloads</h2><ul><li><a href="data:text/csv;base64,$(result["comp_b64"])" download="comparison_table.csv">Download Summary CSV</a></li><li><a href="data:text/csv;base64,$(result["det_b64"])" download="detailed_interactions.csv">Download Detailed CSV</a></li><li><a href="data:image/png;base64,$(result["plot_b64"])" download="residue_interactions.png">Download Chart Image (PNG)</a></li><li><a href="/pymol-script" download="juprot_session.pml">Download PyMOL Session (.pml)</a></li></ul></div><div class=result-section><h2>Residue Interaction Profile (Interaction Type Specific)</h2><img src="data:image/png;base64,$(result["plot_b64"])" alt="Interaction Chart"></div><div class=result-section><h2>Detailed Residue Interactions</h2>$(result["tables_html"])</div><div class=result-section><h2>Analytical Summary</h2><pre>$(ss)</pre></div><div class=footer-link><p><a href=/>Run Another Analysis</a></p></div></div></body></html>"""
#     #     return html(html_body)
#     # end

# route("/analyze", method=POST) do
#         p1 = postpayload(:first_pdb_path); p2 = postpayload(:second_pdb_path)
#         l1 = postpayload(:first_ligand); l2 = postpayload(:second_ligand)
        
#         result = process_pdb_files(p1, p2, l1, l2)
        
#         if haskey(result, "session_data"); SESSION_DATA[] = result["session_data"]; end
#         if haskey(result, "error"); return html("<h1>Error</h1><p>$(result["error"])</p>"); end
        
#         ss = replace(result["summary"], "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;", "\n"=>"<br>")
        
#         html_body = """
#         <!DOCTYPE html><html><head><title>Results</title><link rel="icon" href="/img/favicon.ico"><style>body{font-family:Arial,sans-serif;margin:40px;background:#f4f7f6;color:#333}.container{max-width:900px;margin:auto;background:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.btn{display:inline-block;padding:10px 15px;background:#007bff;color:#fff;border-radius:4px;text-decoration:none;margin-right:10px}.btn:hover{background:#0056b3}h1{text-align:center;color:#0056b3}h2{color:#007bff;border-bottom:1px solid #eee;padding-bottom:5px;margin-top:30px}img{max-width:100%;height:auto;display:block;margin:0 auto}pre{background:#e9ecef;padding:15px;border-radius:5px;white-space:pre-wrap}</style></head><body><div class="container">
#         <div style="text-align:center"><img src="/img/juProt_logo.png" style="height:60px"></div>
#         <h1>Analysis Results</h1>
        
#         <!-- 1. BAR PLOT (Interaction Profile) -->
#         <div>
#             <h2>Interaction Count Profile</h2>
#             <p>Comparison of interaction frequencies per residue across the two complexes.</p>
#             <img src="data:image/png;base64,$(result["plot_b64"])">
#         </div>

#         <!-- 2. ANALYTICAL SUMMARY -->
#         <div>
#             <h2>Analytical Summary</h2>
#             <pre>$ss</pre>
#         </div>

#         <!-- 3. 2D SPATIAL MAP (Real Geometry) -->
#         <div>
#             <h2>2D Spatial Pocket Projection</h2>
#             <p><strong>Geometric Map:</strong> Residues plotted at their <strong>actual X/Y coordinates</strong> relative to the ligand center (Star). This reveals the physical shape of the binding pocket and the spatial distribution of interactions.</p>
#             <img src="data:image/png;base64,$(result["net_b64"])" style="border:1px solid #ddd;padding:10px">
#             <div style="text-align:center;margin-top:10px">
#                 <a href="data:image/png;base64,$(result["net_b64"])" download="spatial_projection.png" class="btn">Download Spatial Map</a>
#             </div>
#         </div>

#         <!-- 4. TABLES -->
#         <div>
#             <h2>Interaction Details</h2>
#             $(result["tables_html"])
#         </div>

#         <!-- 5. DOWNLOADS -->
#         <div style="margin-top:30px;padding-top:20px;border-top:1px solid #eee">
#             <h2>Downloads</h2>
#             <a href="/pymol-script" download="session.pml" class="btn" style="background:#6c757d">Download PyMOL Script</a>
#             <a href="data:text/csv;base64,$(result["comp_b64"])" download="comparison.csv" class="btn">Summary CSV</a>
#             <a href="data:text/csv;base64,$(result["det_b64"])" download="detailed.csv" class="btn">Detailed CSV</a>
#             <a href="data:image/png;base64,$(result["plot_b64"])" download="chart.png" class="btn">Download Bar Chart</a>
#         </div>
        
#         <div style="text-align:center;margin-top:30px"><a href="/">Run Another Analysis</a></div>
#         </div></body></html>
#         """
#         return html(html_body)
#     end
 
#     route("/pymol-script") do
#         if SESSION_DATA[] !== nothing
#             # Generate the script but DO NOT clear the session data
#             return generate_pymol_script(SESSION_DATA[])
#         else
#             # Return a valid PyMOL command that displays an error on screen
#             return "reinitialize; bg_color white; label (0,0,0), \"Error: Session Data Expired. Please Re-run Analysis.\""
#         end
#     end
    
#     route("/how-to-use") do
#      html("""
#      <!DOCTYPE html><html><head><title>juProt: How to Use & Applications</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px 40px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}h1,h2,h3{color:#0056b3}h1{text-align:center;margin-bottom:30px}h2{margin-top:25px;border-bottom:1px solid #eee;padding-bottom:5px}h3{margin-top:20px;color:#007bff}p,li{line-height:1.6}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}ol,ul{padding-left:20px}.footer-link{text-align:center;margin-top:30px}code{background-color:#e9ecef;padding:2px 4px;border-radius:3px;font-family:monospace}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>How to Use juProt & Its Applications</h1><h2>Protocol</h2><ol><li><strong>Prepare PDB Files</strong>: Obtain two protein-ligand complex PDB files (e.g., experimental structures from RCSB PDB or docked complexes from molecular docking experiments). Ensure HETATM records for ligands are present.</li><li><strong>Access the App</strong>: Visit juProt at <code>https://juprot.info/</code>.</li><li><strong>Upload Files</strong>:<ul><li>Upload the first PDB file under "First Complex PDB File".</li><li>Upload the second PDB file under "Second Complex PDB File".</li><li>Click "Load Ligands".</li></ul></li><li><strong>Select Ligands</strong>:<ul><li>juProt will auto-detect potential ligands from each PDB.</li><li>Select the specific ligand of interest from the dropdown menu for each complex.</li><li>Click "Run Analysis".</li></ul></li><li><strong>View Results</strong>:<ul><li>Examine the textual "Analytical Summary" for a quick overview.</li><li>View the "Residue Interaction Plot" for a visual comparison of <strong>interaction frequencies</strong> per residue.</li><li>Download the "Comparison Table (CSV)" for a structured summary of differences and commonalities.</li><li>Download "Detailed Interactions (CSV)" for a complete list of <strong>all interactions</strong> for both complexes with their geometric parameters.</li></ul></li></ol><h2>Applications & Use Cases for juProt</h2><p>juProt is designed to provide rapid comparative insights into protein-ligand interactions, with a <strong>comprehensive focus on the full interactome (Hydrogen bonds, Hydrophobic contacts, Salt Bridges, Pi-Stacking, etc.)</strong>.</p><h3>1. Understanding the Impact of Mutations</h3><p><strong>Scenario:</strong> You have a wild-type protein-ligand structure and a mutant form (e.g., from a SNP or site-directed mutagenesis) bound to the same ligand.</p><p><strong>How juProt Helps:</strong> Compare the native-ligand and mutant-ligand complexes. juProt highlights how the mutation alters the <strong>interaction network</strong>, which can help explain changes in binding affinity, drug efficacy, or resistance mechanisms. (e.g., comparing a wild-type kinase-inhibitor complex with a gatekeeper mutant-inhibitor complex).</p><h3>2. Comparing Different Ligands to the Same Target</h3><p><strong>Scenario:</strong> You have several drug candidates or chemical probes binding to the same protein target.</p><p><strong>How juProt Helps:</strong> Compare Protein+LigandA with Protein+LigandB. juProt helps identify which ligand forms more/different <strong>interactions</strong> and which residues are key common or unique <strong>interaction partners</strong>, aiding in SAR studies and lead optimization.</p><h3>3. Analyzing Ligand Binding to Different Protein Conformations or Isoforms</h3><p><strong>Scenario:</strong> A protein exists in different states (e.g., active/inactive) or as different isoforms, and you have structures of a ligand bound to these variants.</p><p><strong>How juProt Helps:</strong> Compare ProteinState1+Ligand with ProteinState2+Ligand. juProt can reveal how protein structural changes influence the <strong>interaction profile</strong> with a common ligand.</p><h3>4. Validating Docking Poses</h3><p><strong>Scenario:</strong> You have multiple potential binding poses for a ligand from molecular docking.</p><p><strong>How juProt Helps:</strong> Compare the <strong>interaction profile</strong> of different docked poses or a docked pose against an experimental structure (if available) to assess consistency.</p><h3>5. Educational Purposes</h3><p><strong>Scenario:</strong> Teaching students about protein-ligand interactions.</p><p><strong>How juProt Helps:</strong> Provides an easy-to-use tool for students to explore <strong>differential interactions</strong> without complex software or scripting.</p><h3>Future Enhancements</h3><p>juProt is an actively developed open-source project. Future versions aim to include:</p><ul><li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li><li>Allowing analysis of more than two complexes.</li><li>User-configurable parameters for interaction detection.</li></ul><div class=footer-link><p><a href=/>Back to Home</a></p></div></div></body></html>
#      """)
#     end
 
#     route("/about") do
#      html("""
#      <!DOCTYPE html><html><head><title>About juProt</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px 40px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}h1,h2{color:#0056b3}h1{text-align:center;margin-bottom:30px}h2{margin-top:25px;border-bottom:1px solid #eee;padding-bottom:5px}p,li{line-height:1.6}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}.footer-link{text-align:center;margin-top:30px}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>About juProt</h1><p>juProt (<code>https://juprot.info/</code>) is an open-source web application designed to facilitate the rapid and user-friendly comparative analysis of protein-ligand interaction networks, covering the complete spectrum of non-covalent interactions including Hydrogen Bonds, Hydrophobic Contacts, Salt Bridges, Water Bridges, Pi-Stacking, and Halogen Bonds. Understanding how ligands interact with their protein targets, and how these interactions change due to mutations or when comparing different ligands, is fundamental in structural biology, bioinformatics, and drug discovery.</p><h2>Motivation</h2><p>While powerful tools exist for analyzing interactions in a single protein-ligand complex, comparing these interactions across two different complexes often requires manual data extraction, scripting, and collation of results. juProt aims to simplify this process, making comparative interaction analysis accessible to a broader range of researchers and students without requiring extensive computational expertise.</p><h2>Core Technology</h2><p>juProt is developed using the Julia programming language, leveraging the high-performance Genie.jl web framework for its backend and user interface. The core interaction detection is powered by the well-established Protein-Ligand Interaction Profiler (<a href="https://doi.org/10.1093/nar/gkv315" target=_blank>Salentin et al., 2015</a>), a Python-based tool. juProt interfaces with PLIP using the PythonCall.jl package, allowing seamless integration of PLIP's robust algorithms.</p><h2>Current Features</h2><ul><li>Upload of two PDB files for comparison.</li><li>Automated identification of potential ligands with user selection.</li><li>Detection and quantification of full interactomes (H-bonds, Hydrophobic, Salt Bridges, etc.) for each complex.</li><li>Generation of:<ul><li>A comparative summary table (CSV) highlighting common and differential interactions and residues.</li><li>A detailed list of all interactions for both complexes (CSV).</li><li>A bar chart visually comparing interaction profiles per residue (PNG).</li><li>An on-page analytical summary of key findings.</li></ul></li></ul><h2>Future Development</h2><p>juProt is an ongoing project. Future development plans include:</p><ul><li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li><li>Allowing analysis of more than two complexes.</li><li>User-configurable parameters for interaction detection.</li></ul><h2>Development Team</h2><p>juProt was conceived and developed by:<br><a href="https://www.drpaul.cc/" target=_blank>Dr. Benedict Christopher Paul</a><br><a href="https://deepakshankar810.github.io/portfolio/" target=_blank>Deepak S P</a>, MSc Biotechnology<br><a href="https://siva1106.github.io/website/" target=_blank>Siva V</a>, MSc Biotechnology<br>Surya Sekaran, [PhD]</p><p>We also acknowledge the developers of the core libraries used in juProt, including Julia, Genie.jl, PythonCall.jl, PLIP, OpenBabel, and Plots.jl.</p><h2>Open Source & Citation</h2><p>juProt is an open-source project. The source code is available on GitHub at <a href="https://github.com/drbenedictpaul/juprot" target=_blank>https://github.com/drbenedictpaul/juprot</a>.</p><p>We encourage contributions and feedback from the community.</p><p>If you use juProt in your research, please cite:<br><em>[Manuscript is in communication. For now, please cite our GitHub repository.]</em></p><h2>Contact/Feedback</h2><p>For questions, suggestions, or to report issues, please visit our GitHub issues page at <a href="https://github.com/drbenedictpaul/juprot/issues" target=_blank>GitHub Issues</a> or contact <a href="mailto:benedictpaulc@sriramachandra.edu.in">benedictpaulc@sriramachandra.edu.in</a>.</p><div class=footer-link><p><a href=/>Back to Home</a></p></div></div></body></html>
#      """)
#     end
 
#  end # module JuProtGUI



# Filename: app.jl

module JuProtGUI
    # 1. FORCE SILENT MODE (Prevents Cloud Run Crash)
    ENV["GKSwstype"] = "100"

    using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
    using CSV, DataFrames, PythonCall, Printf, Downloads, Plots, Base64

    # --- CONFIGURATION ---
    const LIB_PATH = joinpath(@__DIR__, "..", "lib")
    
    # We use /tmp for temporary generation
    const OUTPUT_DIR = "/tmp/outputs"
    mkpath(OUTPUT_DIR)
    
    include(joinpath(LIB_PATH, "output_results.jl"))
    include(joinpath(LIB_PATH, "detect_hbonds.jl"))
    
    function get_ligand_names(p)
        try
            if !isfile(p); return ["ERROR: File missing at $p"]; end
            plip=pyimport("plip.structure.preparation"); mol=plip.PDBComplex(); mol.load_pdb(p); l=Set{String}()
            if PythonCall.pyhasattr(mol,"ligands")&&mol.ligands!==PythonCall.pybuiltins.None
                for r in mol.ligands
                    if PythonCall.pyhasattr(r,"hetid"); n_py=r.hetid; if n_py!==PythonCall.pybuiltins.None; n=pyconvert(String,n_py); if !isempty(n)&&n!="HOH"; push!(l,n); end; end; end
                end
            end
            if isempty(l)&&PythonCall.pyhasattr(mol,"res_het")
                for r in mol.res_het; if PythonCall.pyhasattr(r,"hetid"); n=pyconvert(String,r.hetid); if !isempty(n)&&n!="HOH"; push!(l,n); end; end; end
            end
            if isempty(l); return ["No Ligands Detected"]; end
            return collect(l)
        catch e; err_msg = sprint(showerror, e); return ["ERROR: $err_msg"]; end
    end
    
    function process_pdb_files(p1, p2, l1, l2)
        cr = Dict{String, Dict}()
        fk = "first_complex.pdb"
        sk = "second_complex.pdb"
        files = [(path=p1, ligand=l1, key=fk), (path=p2, ligand=l2, key=sk)]
        ok = true
        
        for i in files
            try
                (ints, lig_meta) = detect_hbonds(i.path, i.ligand)
                ir = Dict{String, Int}()
                for d in ints
                     c_display = isempty(d.prot_chain) ? "_" : d.prot_chain
                     rk = "$(c_display) $(d.prot_resn) $(d.prot_resi)"
                     ir[rk] = get(ir, rk, 0) + 1
                end
                cr[i.key] = Dict(
                    :hbonds => ints,
                    :interacting_residues => ir,
                    :ligand_resname => i.ligand,
                    :ligand_data => lig_meta,
                    :original_path => i.path
                )
            catch e; println("Error processing $(i.path): $e"); ok=false; end
        end
        
        if !ok; return Dict("error" => "Processing failed."); end
        
        session_id = string(abs(rand(Int)))[1:6]
        ocsv = joinpath(OUTPUT_DIR, "comparison_table_$(session_id).csv")
        dcsv = joinpath(OUTPUT_DIR, "detailed_interactions_$(session_id).csv")
        ppath = joinpath(OUTPUT_DIR, "residue_interactions_$(session_id).png")
        net_path = joinpath(OUTPUT_DIR, "network_map_$(session_id).png")
 
        try
            generate_stacked_chart(cr, fk, sk, ppath)
            generate_2d_network_map(cr, fk, sk, net_path)
            compare_and_save_results(cr, fk, sk, ocsv, dcsv, "dummy.png") 
            st = capture_analytical_summary(cr, fk, sk)
            tables_html = generate_interaction_tables(cr, fk, sk)
            
            comp_b64 = base64encode(read(ocsv))
            det_b64 = base64encode(read(dcsv))
            plot_b64 = base64encode(read(ppath))
            net_b64 = base64encode(read(net_path))
            
            rm(ocsv, force=true); rm(dcsv, force=true); rm(ppath, force=true); rm(net_path, force=true)
 
            return Dict(
                "comp_b64" => comp_b64, 
                "det_b64" => det_b64, 
                "plot_b64" => plot_b64, 
                "net_b64" => net_b64, 
                "summary" => st, 
                "tables_html" => tables_html,
                "session_data" => cr # This might be removable if nothing else uses it.
            )
        catch e; return Dict("error" => "Visualization failed: $(sprint(showerror, e))"); end
    end
    
    function capture_analytical_summary(cr,fk,sk);try;io=IOBuffer();print_analytical_summary(io,cr,fk,sk);return String(take!(io));catch;return"Summary failed.";end;end
 
    route("/") do
     html("""
     <!DOCTYPE html><html><head><title>juProt</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon">
     <style>
        body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}
        .container{max-width:800px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}
        .form-group{margin-bottom:20px}
        label{display:block;margin-bottom:8px;font-weight:700}
        input[type=file],input[type=text]{width:calc(100% - 22px);padding:10px;border:1px solid #ccc;border-radius:4px}
        button{padding:10px 20px;background:#007bff;color:#fff;border:none;cursor:pointer;border-radius:4px;font-size:16px}
        button:hover{background:#0056b3}
        h1{color:#0056b3;text-align:center}
        a{color:#007bff;text-decoration:none}
        .footer-links-container{display:flex;justify-content:center;gap:30px;margin-top:40px;padding-top:20px;border-top:1px solid #eee}
        .divider{text-align:center;font-weight:700;color:#aaa;margin:20px 0}
        .validation-popup{display:none;position:absolute;top:35px;left:110px;background-color:#fff;border:1px solid #ccc;border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,.15);padding:8px 12px;z-index:100;white-space:nowrap;font-size:14px;color:#333}.validation-popup::before{content:'';position:absolute;bottom:100%;left:20px;border-width:7px;border-style:solid;border-color:transparent transparent #ccc transparent}.validation-popup::after{content:'';position:absolute;bottom:100%;left:21px;border-width:6px;border-style:solid;border-color:transparent transparent #fff transparent}.validation-popup .icon{display:inline-block;background-color:#ff8552;color:#fff;width:16px;height:16px;border-radius:3px;text-align:center;font-weight:700;line-height:16px;margin-right:8px;font-size:12px}
        
        /* LOADING OVERLAY CSS */
        #loader {
            display: none;
            position: fixed; top: 0; left: 0; width: 100%; height: 100%;
            background: rgba(255,255,255,0.9);
            z-index: 1000;
            flex-direction: column; justify-content: center; align-items: center;
        }
        .spinner {
            border: 8px solid #f3f3f3; border-top: 8px solid #007bff;
            border-radius: 50%; width: 60px; height: 60px;
            animation: spin 1s linear infinite;
        }
        @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
        .loading-text { margin-top: 20px; font-size: 18px; color: #0056b3; font-weight: bold; }
     </style>
     <script>
        function validateForm(){
            document.getElementById('alert_1').style.display='none';
            document.getElementById('alert_2').style.display='none';
            var f1=document.getElementById('first_pdb').value;
            var i1=document.getElementById('first_pdb_id').value.trim();
            var c1=(f1!==""||i1!=="");
            var f2=document.getElementById('second_pdb').value;
            var i2=document.getElementById('second_pdb_id').value.trim();
            var c2=(f2!==""||i2!=="");
            var v=true;
            if(!c1){document.getElementById('alert_1').style.display='flex';v=false}
            if(!c2){document.getElementById('alert_2').style.display='flex';v=false}
            
            // SHOW LOADER IF VALID
            if(v) {
                document.getElementById('loader').style.display = 'flex';
            }
            return v;
        }
     </script>
     </head>
     <body>
        <!-- LOADER DIV -->
        <div id="loader">
            <div class="spinner"></div>
            <div class="loading-text">Uploading and Parsing Structures...<br><span style="font-size:14px;color:#666;font-weight:normal">This may take a few seconds.</span></div>
        </div>

        <div class=container>
            <div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div>
            <h1>juProt: Protein-Ligand Interaction Analyzer</h1>
            <p style="text-align:center;color:#666">Comparative analysis of the complete protein-ligand interactome.</p>
            <form action=/select-ligands method=post enctype=multipart/form-data onsubmit="return validateForm()">
                <div class=form-group>
                    <label for=first_pdb>Upload First Complex (PDB File):</label>
                    <input type=file id=first_pdb name=first_pdb accept=.pdb>
                    <div id=alert_1 class=validation-popup><span class=icon>!</span>Please select a file or enter a PDB ID.</div>
                </div>
                <div class=form-group>
                    <label for=second_pdb>Upload Second Complex (PDB File):</label>
                    <input type=file id=second_pdb name=second_pdb accept=.pdb>
                    <div id=alert_2 class=validation-popup><span class=icon>!</span>Please select a file or enter a PDB ID.</div>
                </div>
                <div class=divider>OR</div>
                <div class=form-group><label for=first_pdb_id>Enter First PDB ID (e.g., 3EQM):</label><input type=text id=first_pdb_id name=first_pdb_id placeholder="PDB ID"></div>
                <div class=form-group><label for=second_pdb_id>Enter Second PDB ID:</label><input type=text id=second_pdb_id name=second_pdb_id placeholder="PDB ID"></div>
                <button type=submit>Load Ligands</button>
            </form>
            <div class=footer-links-container>
                <p><a href=/how-to-use>How to Use & Applications</a></p>
                <p><a href=/about>About juProt</a></p>
            </div>
        </div>
     </body></html>
     """)
    end
    
    route("/select-ligands", method=POST) do
        first_payload=haskey(filespayload(),"first_pdb") ? filespayload("first_pdb") : nothing;
        second_payload=haskey(filespayload(),"second_pdb") ? filespayload("second_pdb") : nothing
        f_id=postpayload(:first_pdb_id,""); s_id=postpayload(:second_pdb_id,"")
        mkpath(OUTPUT_DIR); tp1=""; tp2=""; fn1=""; fn2=""
        try
            if !isempty(f_id); tp1=joinpath(OUTPUT_DIR,"$(f_id).pdb"); Downloads.download("https://files.rcsb.org/download/$(f_id).pdb",tp1); fn1="$(f_id).pdb"
            elseif first_payload!==nothing&&!isempty(first_payload.data); tp1=joinpath(OUTPUT_DIR,"temp_first_"*string(randn())[3:8]*".pdb"); write(tp1,first_payload.data); fn1=first_payload.name
            else; return html("<h1>Error</h1><p>No input for first complex.</p>"); end
            if !isempty(s_id); tp2=joinpath(OUTPUT_DIR,"$(s_id).pdb"); Downloads.download("https://files.rcsb.org/download/$(s_id).pdb",tp2); fn2="$(s_id).pdb"
            elseif second_payload!==nothing&&!isempty(second_payload.data); tp2=joinpath(OUTPUT_DIR,"temp_second_"*string(randn())[3:8]*".pdb"); write(tp2,second_payload.data); fn2=second_payload.name
            else; return html("<h1>Error</h1><p>No input for second complex.</p>"); end
            lig1=get_ligand_names(tp1); lig2=get_ligand_names(tp2)
            if isempty(lig1)&&isempty(lig2); return html("<h1>Error</h1><p>No ligands found.</p>"); end
            opt1=join(["<option value='$l'>$l</option>" for l in lig1],""); opt2=join(["<option value='$l'>$l</option>" for l in lig2],"")
            
            html_body="""<!DOCTYPE html><html><head><title>juProt: Select Ligands</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon">
            <style>
                body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}
                .container{max-width:800px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}
                .form-group{margin-bottom:20px}label{display:block;margin-bottom:8px;font-weight:700}
                select{width:100%;padding:10px;border:1px solid #ccc;border-radius:4px;box-sizing:border-box}
                button{padding:10px 20px;background:#007bff;color:#fff;border:none;cursor:pointer;border-radius:4px;font-size:16px}
                button:hover{background:#0056b3}h1{color:#0056b3;text-align:center}.footer-link{text-align:center;margin-top:20px}a{color:#007bff;text-decoration:none}
                #loader {
                    display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%;
                    background: rgba(255,255,255,0.9); z-index: 1000;
                    flex-direction: column; justify-content: center; align-items: center;
                }
                .spinner {
                    border: 8px solid #f3f3f3; border-top: 8px solid #007bff;
                    border-radius: 50%; width: 60px; height: 60px;
                    animation: spin 1s linear infinite;
                }
                @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
                .loading-text { margin-top: 20px; font-size: 18px; color: #0056b3; font-weight: bold; }
            </style>
            <script> function showLoader() { document.getElementById('loader').style.display = 'flex'; } </script>
            </head><body>
                <div id="loader">
                    <div class="spinner"></div>
                    <div class="loading-text">Calculating Interactions...<br><span style="font-size:14px;color:#666;font-weight:normal">This involves complex calculations (H-bonds, Pi-Stacking, etc).<br>Please wait...</span></div>
                </div>
                <div class=container>
                    <div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div>
                    <h1>Select Target Ligands</h1>
                    <form action=/analyze method=post onsubmit="showLoader()">
                        <input type=hidden name=first_pdb_path value="$tp1">
                        <input type=hidden name=second_pdb_path value="$tp2">
                        <div class=form-group>
                            <label>First Complex Ligand ($fn1):</label>
                            <select name=first_ligand required>$opt1</select>
                        </div>
                        <div class=form-group>
                            <label>Second Complex Ligand ($fn2):</label>
                            <select name=second_ligand required>$opt2</select>
                        </div>
                        <button type=submit>Run Analysis</button>
                    </form>
                    <div class=footer-link><p><a href=/>Back</a></p></div>
                </div>
            </body></html>"""
            return html(html_body)
        catch e; return html("<h1>Error</h1><p>$(sprint(showerror,e))</p>"); end
    end

    route("/analyze", method=POST) do
        p1 = postpayload(:first_pdb_path); p2 = postpayload(:second_pdb_path)
        l1 = postpayload(:first_ligand); l2 = postpayload(:second_ligand)
        
        result = process_pdb_files(p1, p2, l1, l2)
        
        if haskey(result, "error"); return html("<h1>Error</h1><p>$(result["error"])</p>"); end
        
        ss = replace(result["summary"], "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;", "\n"=>"<br>")
        
        html_body = """
        <!DOCTYPE html><html><head><title>Results</title><link rel="icon" href="/img/favicon.ico"><style>body{font-family:Arial,sans-serif;margin:40px;background:#f4f7f6;color:#333}.container{max-width:900px;margin:auto;background:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.btn{display:inline-block;padding:10px 15px;background:#007bff;color:#fff;border-radius:4px;text-decoration:none;margin-right:10px}.btn:hover{background:#0056b3}h1{text-align:center;color:#0056b3}h2{color:#007bff;border-bottom:1px solid #eee;padding-bottom:5px;margin-top:30px}img{max-width:100%;height:auto;display:block;margin:0 auto}pre{background:#e9ecef;padding:15px;border-radius:5px;white-space:pre-wrap}</style></head><body><div class="container">
        <div style="text-align:center"><img src="/img/juProt_logo.png" style="height:60px"></div>
        <h1>Analysis Results</h1>
        <div><h2>Interaction Count Profile</h2><p>Comparison of interaction frequencies per residue across the two complexes.</p><img src="data:image/png;base64,$(result["plot_b64"])"></div>
        <div><h2>Analytical Summary</h2><pre>$ss</pre></div>
        <div><h2>2D Spatial Pocket Projection</h2><p><strong>Geometric Map:</strong> Residues plotted at their <strong>actual X/Y coordinates</strong> relative to the ligand center (Star). This reveals the physical shape of the binding pocket and the spatial distribution of interactions.</p><img src="data:image/png;base64,$(result["net_b64"])" style="border:1px solid #ddd;padding:10px"><div style="text-align:center;margin-top:10px"><a href="data:image/png;base64,$(result["net_b64"])" download="spatial_projection.png" class="btn">Download Spatial Map</a></div></div>
        <div><h2>Interaction Details</h2>$(result["tables_html"])</div>
        <div style="margin-top:30px;padding-top:20px;border-top:1px solid #eee">
            <h2>Downloads</h2>
            <a href="data:text/csv;base64,$(result["comp_b64"])" download="comparison.csv" class="btn">Summary CSV</a>
            <a href="data:text/csv;base64,$(result["det_b64"])" download="detailed.csv" class="btn">Detailed CSV</a>
            <a href="data:image/png;base64,$(result["plot_b64"])" download="chart.png" class="btn">Download Bar Chart</a>
        </div>
        <div style="text-align:center;margin-top:30px"><a href="/">Run Another Analysis</a></div>
        </div></body></html>
        """
        return html(html_body)
    end
    
    route("/how-to-use") do
     html("""
     <!DOCTYPE html><html><head><title>juProt: How to Use & Applications</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px 40px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}h1,h2,h3{color:#0056b3}h1{text-align:center;margin-bottom:30px}h2{margin-top:25px;border-bottom:1px solid #eee;padding-bottom:5px}h3{margin-top:20px;color:#007bff}p,li{line-height:1.6}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}ol,ul{padding-left:20px}.footer-link{text-align:center;margin-top:30px}code{background-color:#e9ecef;padding:2px 4px;border-radius:3px;font-family:monospace}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>How to Use juProt & Its Applications</h1><h2>Protocol</h2><ol><li><strong>Prepare PDB Files</strong>: Obtain two protein-ligand complex PDB files (e.g., experimental structures from RCSB PDB or docked complexes from molecular docking experiments). Ensure HETATM records for ligands are present.</li><li><strong>Access the App</strong>: Visit juProt at <code>https://juprot.info/</code>.</li><li><strong>Upload Files</strong>:<ul><li>Upload the first PDB file under "First Complex PDB File".</li><li>Upload the second PDB file under "Second Complex PDB File".</li><li>Click "Load Ligands".</li></ul></li><li><strong>Select Ligands</strong>:<ul><li>juProt will auto-detect potential ligands from each PDB.</li><li>Select the specific ligand of interest from the dropdown menu for each complex.</li><li>Click "Run Analysis".</li></ul></li><li><strong>View Results</strong>:<ul><li>Examine the textual "Analytical Summary" for a quick overview.</li><li>View the "Residue Interaction Plot" for a visual comparison of <strong>interaction frequencies</strong> per residue.</li><li>Download the "Comparison Table (CSV)" for a structured summary of differences and commonalities.</li><li>Download "Detailed Interactions (CSV)" for a complete list of <strong>all interactions</strong> for both complexes with their geometric parameters.</li></ul></li></ol><h2>Applications & Use Cases for juProt</h2><p>juProt is designed to provide rapid comparative insights into protein-ligand interactions, with a <strong>comprehensive focus on the full interactome (Hydrogen bonds, Hydrophobic contacts, Salt Bridges, Pi-Stacking, etc.)</strong>.</p><h3>1. Understanding the Impact of Mutations</h3><p><strong>Scenario:</strong> You have a wild-type protein-ligand structure and a mutant form (e.g., from a SNP or site-directed mutagenesis) bound to the same ligand.</p><p><strong>How juProt Helps:</strong> Compare the native-ligand and mutant-ligand complexes. juProt highlights how the mutation alters the <strong>interaction network</strong>, which can help explain changes in binding affinity, drug efficacy, or resistance mechanisms. (e.g., comparing a wild-type kinase-inhibitor complex with a gatekeeper mutant-inhibitor complex).</p><h3>2. Comparing Different Ligands to the Same Target</h3><p><strong>Scenario:</strong> You have several drug candidates or chemical probes binding to the same protein target.</p><p><strong>How juProt Helps:</strong> Compare Protein+LigandA with Protein+LigandB. juProt helps identify which ligand forms more/different <strong>interactions</strong> and which residues are key common or unique <strong>interaction partners</strong>, aiding in SAR studies and lead optimization.</p><h3>3. Analyzing Ligand Binding to Different Protein Conformations or Isoforms</h3><p><strong>Scenario:</strong> A protein exists in different states (e.g., active/inactive) or as different isoforms, and you have structures of a ligand bound to these variants.</p><p><strong>How juProt Helps:</strong> Compare ProteinState1+Ligand with ProteinState2+Ligand. juProt can reveal how protein structural changes influence the <strong>interaction profile</strong> with a common ligand.</p><h3>4. Validating Docking Poses</h3><p><strong>Scenario:</strong> You have multiple potential binding poses for a ligand from molecular docking.</p><p><strong>How juProt Helps:</strong> Compare the <strong>interaction profile</strong> of different docked poses or a docked pose against an experimental structure (if available) to assess consistency.</p><h3>5. Educational Purposes</h3><p><strong>Scenario:</strong> Teaching students about protein-ligand interactions.</p><p><strong>How juProt Helps:</strong> Provides an easy-to-use tool for students to explore <strong>differential interactions</strong> without complex software or scripting.</p><h3>Future Enhancements</h3><p>juProt is an actively developed open-source project. Future versions aim to include:</p><ul><li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li><li>Allowing analysis of more than two complexes.</li><li>User-configurable parameters for interaction detection.</li></ul><div class=footer-link><p><a href=/>Back to Home</a></p></div></div></body></html>
     """)
    end
 
    route("/about") do
     html("""
     <!DOCTYPE html><html><head><title>About juProt</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px 40px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}h1,h2{color:#0056b3}h1{text-align:center;margin-bottom:30px}h2{margin-top:25px;border-bottom:1px solid #eee;padding-bottom:5px}p,li{line-height:1.6}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}.footer-link{text-align:center;margin-top:30px}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>About juProt</h1><p>juProt (<code>https://juprot.info/</code>) is an open-source web application designed to facilitate the rapid and user-friendly comparative analysis of protein-ligand interaction networks, covering the complete spectrum of non-covalent interactions including Hydrogen Bonds, Hydrophobic Contacts, Salt Bridges, Water Bridges, Pi-Stacking, and Halogen Bonds. Understanding how ligands interact with their protein targets, and how these interactions change due to mutations or when comparing different ligands, is fundamental in structural biology, bioinformatics, and drug discovery.</p><h2>Motivation</h2><p>While powerful tools exist for analyzing interactions in a single protein-ligand complex, comparing these interactions across two different complexes often requires manual data extraction, scripting, and collation of results. juProt aims to simplify this process, making comparative interaction analysis accessible to a broader range of researchers and students without requiring extensive computational expertise.</p><h2>Core Technology</h2><p>juProt is developed using the Julia programming language, leveraging the high-performance Genie.jl web framework for its backend and user interface. The core interaction detection is powered by the well-established Protein-Ligand Interaction Profiler (<a href="https://doi.org/10.1093/nar/gkv315" target=_blank>Salentin et al., 2015</a>), a Python-based tool. juProt interfaces with PLIP using the PythonCall.jl package, allowing seamless integration of PLIP's robust algorithms.</p><h2>Current Features</h2><ul><li>Upload of two PDB files for comparison.</li><li>Automated identification of potential ligands with user selection.</li><li>Detection and quantification of full interactomes (H-bonds, Hydrophobic, Salt Bridges, etc.) for each complex.</li><li>Generation of:<ul><li>A comparative summary table (CSV) highlighting common and differential interactions and residues.</li><li>A detailed list of all interactions for both complexes (CSV).</li><li>A bar chart visually comparing interaction profiles per residue (PNG).</li><li>An on-page analytical summary of key findings.</li></ul></li></ul><h2>Future Development</h2><p>juProt is an ongoing project. Future development plans include:</p><ul><li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li><li>Allowing analysis of more than two complexes.</li><li>User-configurable parameters for interaction detection.</li></ul><h2>Development Team</h2><p>juProt was conceived and developed by:<br><a href="https://www.drpaul.cc/" target=_blank>Dr. Benedict Christopher Paul</a><br><a href="https://deepakshankar810.github.io/portfolio/" target=_blank>Deepak S P</a>, MSc Biotechnology<br><a href="https://siva1106.github.io/website/" target=_blank>Siva V</a>, MSc Biotechnology<br>Surya Sekaran, [PhD]</p><p>We also acknowledge the developers of the core libraries used in juProt, including Julia, Genie.jl, PythonCall.jl, PLIP, OpenBabel, and Plots.jl.</p><h2>Open Source & Citation</h2><p>juProt is an open-source project. The source code is available on GitHub at <a href="https://github.com/drbenedictpaul/juprot" target=_blank>https://github.com/drbenedictpaul/juprot</a>.</p><p>We encourage contributions and feedback from the community.</p><p>If you use juProt in your research, please cite:<br><em>[Manuscript is in communication. For now, please cite our GitHub repository.]</em></p><h2>Contact/Feedback</h2><p>For questions, suggestions, or to report issues, please visit our GitHub issues page at <a href="https://github.com/drbenedictpaul/juprot/issues" target=_blank>GitHub Issues</a> or contact <a href="mailto:benedictpaulc@sriramachandra.edu.in">benedictpaulc@sriramachandra.edu.in</a>.</p><div class=footer-link><p><a href=/>Back to Home</a></p></div></div></body></html>
     """)
    end
 
 end # module JuProtGUI