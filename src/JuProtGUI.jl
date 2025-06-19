# JuProtGUI.jl (located in your src/ folder)

module JuProtGUI
   using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
   using CSV, DataFrames, Plots, PythonCall

   # --- Debugging (can be removed once includes are confirmed working) ---
   println("Current file's directory (__DIR__): ", @__DIR__)
   path_to_output_results = joinpath(@__DIR__, "..", "lib", "output_results.jl")
   path_to_detect_hbonds = joinpath(@__DIR__, "..", "lib", "detect_hbonds.jl")
   println("Attempting to include output_results.jl from: ", path_to_output_results)
   println("File exists? ", isfile(path_to_output_results))
   println("Attempting to include detect_hbonds.jl from: ", path_to_detect_hbonds)
   println("File exists? ", isfile(path_to_detect_hbonds))
   # --- End Debugging ---

   include(path_to_output_results)
   include(path_to_detect_hbonds)

   const OUTPUT_DIR = "public/outputs"

   function get_ligand_names(pdb_file_path::String)
       println("--- Starting get_ligand_names for: '$pdb_file_path' ---")
       try
           plip = pyimport("plip.structure.preparation")
           mol = plip.PDBComplex()
           mol.load_pdb(pdb_file_path)
           println("DEBUG get_ligand_names: PDB loaded: '$pdb_file_path'")
           mol.analyze()
           println("DEBUG get_ligand_names: mol.analyze() called.")

           ligands = Set{String}()
           if mol.ligands !== PythonCall.pybuiltins.None && pylen(mol.ligands) > 0
               println("DEBUG get_ligand_names: mol.ligands is not None and not empty. Length: $(pylen(mol.ligands))")
               for res_obj in mol.ligands
                   println("DEBUG get_ligand_names: Processing res_obj: $(pystr(res_obj))")
                   if res_obj !== PythonCall.pybuiltins.None && PythonCall.pyhasattr(res_obj, "hetid")
                       lig_name_py = res_obj.hetid
                       println("DEBUG get_ligand_names: res_obj.hetid (raw PyObject): $(pystr(lig_name_py))")
                       if lig_name_py !== PythonCall.pybuiltins.None
                           lig_name = pyconvert(String, lig_name_py)
                           println("DEBUG get_ligand_names: Converted hetid to Julia String: '$lig_name'")
                           if !isempty(lig_name) && lig_name != "HOH"
                               push!(ligands, lig_name)
                               println("DEBUG get_ligand_names: Added '$lig_name' to ligands set.")
                           else
                               println("DEBUG get_ligand_names: Ligand '$lig_name' is empty or HOH. Skipping.")
                           end
                       else
                           println("DEBUG get_ligand_names: res_obj.hetid was Python None. Skipping.")
                       end
                   else
                       obj_str = res_obj === PythonCall.pybuiltins.None ? "Python None" : pystr(res_obj)
                       if res_obj === PythonCall.pybuiltins.None
                           println("DEBUG get_ligand_names: res_obj itself is Python None. Skipping.")
                       elseif !PythonCall.pyhasattr(res_obj, "hetid") 
                           println("DEBUG get_ligand_names: res_obj ($(obj_str)) does not have attribute 'hetid'. Skipping.")
                       else # Should not happen if the outer if is structured well
                           println("DEBUG get_ligand_names: Unhandled case for res_obj ($(obj_str)). Skipping.")
                       end
                   end
               end
           else
               list_len_str = mol.ligands === PythonCall.pybuiltins.None ? "Python None" : string(pylen(mol.ligands))
               println("DEBUG get_ligand_names: mol.ligands is $(list_len_str == "Python None" ? "None" : "empty (length $list_len_str)").")
           end
           
           if isempty(ligands)
               println("FINAL get_ligand_names: No ligands collected for '$pdb_file_path'. Returning empty list.")
           else
               println("FINAL get_ligand_names: Detected potential ligands in '$pdb_file_path': ", join(ligands, ", "))
           end
           return collect(ligands)
       catch e
           println("Error in get_ligand_names for PDB file $pdb_file_path:")
           Base.showerror(stdout, e)
           Base.show_backtrace(stdout, Base.catch_backtrace())
           println()

           if e isa PythonCall.PyException
               println("--- Wrapped Python Exception in get_ligand_names ---")
               py_exc_obj = e.exc
               py_traceback_mod = pyimport("traceback")
               println("Python Exception (pyrepr): $(PythonCall.pyrepr(String, py_exc_obj))")
               println("Attempting to print Python traceback using traceback.print_exception():")
               try
                   py_traceback_mod.print_exception(pytype(py_exc_obj), py_exc_obj, py_exc_obj.__traceback__)
               catch py_print_err
                   println("Failed to print Python traceback using traceback.print_exception(): $py_print_err")
                   try
                       formatted_tb = py_traceback_mod.format_exc()
                       println("Fallback format_exc():\n", pyconvert(String, formatted_tb))
                   catch py_format_err
                       println("format_exc() also failed: $py_format_err")
                   end
               end
               println("-------------------------------------------------")
           else
                println("The caught exception in get_ligand_names was a Julia exception, not directly a Python one.")
           end
           return String[]
       end
   end

   function process_pdb_files(first_pdb_path::String, second_pdb_path::String, first_ligand::String, second_ligand::String)
       println("Starting process_pdb_files with files: ", first_pdb_path, ", ", second_pdb_path)
       println("Ligands selected: First=", first_ligand, ", Second=", second_ligand)
       complex_results = Dict{String, Dict}()
       
       mkpath(OUTPUT_DIR)
       
       first_complex_key = "first_complex.pdb"
       second_complex_key = "second_complex.pdb"

       pdb_files_to_process = [(path=first_pdb_path, ligand=first_ligand, key=first_complex_key),
                               (path=second_pdb_path, ligand=second_ligand, key=second_complex_key)]

       all_processed_successfully = true
       for item in pdb_files_to_process
           pdb_file = item.path
           ligand_resname = item.ligand
           complex_key_name = item.key
           try
               println("Processing PDB file: ", pdb_file, " for ligand: ", ligand_resname)
               hbonds = detect_hbonds(pdb_file, ligand_resname) # detect_hbonds returns [] on Python error
               
               # If detect_hbonds itself had a critical Julia error and rethrew, this catch block would get it.
               # Assuming detect_hbonds handles its Python errors and returns [], processing continues.
               
               println("Found ", length(hbonds), " H-bonds for ", ligand_resname, " in ", pdb_file)
               
               interacting_residues = Dict{String, Int}()
               for bond_detail in hbonds
                   res_key = "$(bond_detail.prot_resn) $(bond_detail.prot_resi)"
                   interacting_residues[res_key] = get(interacting_residues, res_key, 0) + 1
               end
               
               max_residue_str = ""
               max_count = 0
               if !isempty(interacting_residues)
                   max_res_tuple = findmax(interacting_residues)
                   max_count = max_res_tuple[1]
                   max_residue_str = max_res_tuple[2]
               end

               complex_results[complex_key_name] = Dict(
                   :hbonds => hbonds,
                   :interacting_residues => interacting_residues,
                   :max_residue => max_residue_str,
                   :max_count => max_count,
                   :ligand_resname => ligand_resname,
                   :original_path => pdb_file
               )
           catch e # Catch Julia errors in this loop's logic
               err_msg = "JULIA ERROR during processing of $pdb_file for $ligand_resname: $(sprint(showerror, e))"
               println(err_msg)
               Base.show_backtrace(stdout, Base.catch_backtrace())
               println()
               complex_results[complex_key_name] = Dict("error_message" => err_msg, :hbonds => [], :interacting_residues => Dict(), :ligand_resname => ligand_resname) # Store error
               all_processed_successfully = false # Mark that at least one PDB failed
               # Continue to next PDB if possible, rather than returning immediately
           end
       end
       
       # If not all PDBs were processed successfully (e.g., Julia error in the loop), return an error.
       # Note: If detect_hbonds had a Python error, it returns [], so it counts as "processed" with 0 bonds here.
       # This `if` is mainly for Julia errors in the loop above.
       if !all_processed_successfully || !haskey(complex_results, first_complex_key) || !haskey(complex_results, second_complex_key)
           # Consolidate error messages if available
           errors = []
           for key_val in [first_complex_key, second_complex_key]
                if haskey(complex_results, key_val) && haskey(complex_results[key_val], "error_message")
                    push!(errors, complex_results[key_val]["error_message"])
                elseif !haskey(complex_results, key_val)
                    push!(errors, "Processing data for $key_val is missing.")
                end
           end
           final_error_msg = isempty(errors) ? "One or both PDB files could not be processed fully (unknown reason)." : join(errors, "; ")
           println("Error: " * final_error_msg)
           return Dict("error" => final_error_msg)
       end

       output_comparison_csv = joinpath(OUTPUT_DIR, "comparison_table.csv")
       plot_file_path = joinpath(OUTPUT_DIR, "residue_interactions.png")
       detailed_interactions_csv = joinpath(OUTPUT_DIR, "detailed_interactions.csv")

       try
           println("Generating comparison and summary using keys: ", keys(complex_results))
           compare_and_save_results(complex_results, first_complex_key, second_complex_key, output_comparison_csv, detailed_interactions_csv, plot_file_path)
           
           summary_text = capture_analytical_summary(complex_results, first_complex_key, second_complex_key)
           println("Captured summary: ", summary_text)
           
           return Dict(
               "comparison_table" => "/outputs/comparison_table.csv",
               "detailed_interactions" => "/outputs/detailed_interactions.csv",
               "plot_file" => "/outputs/residue_interactions.png",
               "summary" => summary_text
           )
       catch e
           err_msg = "Analysis or summary generation failed: $(sprint(showerror, e))"
           println(err_msg)
           Base.show_backtrace(stdout, Base.catch_backtrace())
           println()
           return Dict("error" => err_msg)
       end
   end

   function capture_analytical_summary(complex_results, first_key, second_key)
       println("Capturing analytical summary for keys: ", first_key, ", ", second_key)
       summary_str = ""
       try
           io = IOBuffer()
           print_analytical_summary(io, complex_results, first_key, second_key)
           summary_str = String(take!(io))
           
           if isempty(summary_str)
               println("Warning: Empty analytical summary captured")
               summary_str = "No summary could be generated."
           else
               println("Summary captured successfully.")
           end
           return summary_str
       catch e
           err_msg = "Error capturing summary: $(sprint(showerror, e))"
           println(err_msg)
           return err_msg
       end
   end

   # In src/JuProtGUI.jl
route("/") do
    html("""
    <!DOCTYPE html>
    <html>
    <head>
        <title>juProt</title>
        <link rel="icon" href="/img/favicon.ico" type="image/x-icon">
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; background-color: #f4f7f6; color: #333; }
            .container { max-width: 800px; margin: auto; background-color: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
            .form-group { margin-bottom: 20px; }
            label { display: block; margin-bottom: 8px; font-weight: bold; }
            input[type=file], select { width: calc(100% - 18px); padding: 10px; border: 1px solid #ccc; border-radius: 4px; }
            button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; border-radius: 4px; font-size: 16px; }
            button:hover { background: #0056b3; }
            h1 { color: #0056b3; text-align: center; }
            a { color: #007bff; text-decoration: none; }
            a:hover { text-decoration: underline; }
            .footer-links-container { 
                display: flex; 
                justify-content: center; 
                gap: 30px; /* Increased gap for better separation */
                margin-top: 40px; 
                padding-top: 20px;
                border-top: 1px solid #eee; 
            }
            .footer-links-container p {
                margin: 0; 
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
            <h1>juProt: Protein-Ligand Interaction Analysis</h1>
            <form action="/select-ligands" method="post" enctype="multipart/form-data">
                <div class="form-group">
                    <label for="first_pdb">First Complex PDB File:</label>
                    <input type="file" id="first_pdb" name="first_pdb" accept=".pdb" required>
                </div>
                <div class="form-group">
                    <label for="second_pdb">Second Complex PDB File:</label>
                    <input type="file" id="second_pdb" name="second_pdb" accept=".pdb" required>
                </div>
                <button type="submit">Load Ligands</button>
            </form>

            <div class="footer-links-container">
                <p><a href="/how-to-use">How to Use & Applications</a></p>
                <p><a href="/about">About juProt</a></p>
            </div>

        </div>
    </body>
    </html>
    """)
end

   route("/select-ligands", method=POST) do
       first_file_payload = filespayload("first_pdb")
       second_file_payload = filespayload("second_pdb")

       if first_file_payload === nothing || isempty(first_file_payload.data) ||
          second_file_payload === nothing || isempty(second_file_payload.data)
           return html("<h1>Error</h1><p>Both PDB files must be uploaded.</p><p><a href='/'>Try again</a></p>")
       end
       
       mkpath(OUTPUT_DIR)

       temp_first_path = joinpath(OUTPUT_DIR, "temp_first_" * string(randn())[3:8] * ".pdb")
       temp_second_path = joinpath(OUTPUT_DIR, "temp_second_" * string(randn())[3:8] * ".pdb")

       try
           write(temp_first_path, first_file_payload.data)
           write(temp_second_path, second_file_payload.data)

           first_ligands = get_ligand_names(temp_first_path)
           second_ligands = get_ligand_names(temp_second_path)

           # Cleanup temp files after get_ligand_names, regardless of outcome, unless an error below prevents it
           # If get_ligand_names itself errors hard, this might not be reached.

           if isempty(first_ligands) && isempty(second_ligands)
               msg = "No ligands found in either PDB file. Please check PDB HETATM records or console logs for errors from get_ligand_names."
               println("Error: " * msg)
               isfile(temp_first_path) && rm(temp_first_path, force=true); isfile(temp_second_path) && rm(temp_second_path, force=true)
               return html("<h1>Error</h1><p>$msg</p><p><a href='/'>Try again</a></p>")
           elseif isempty(first_ligands)
               msg = "No ligands found in the first PDB file ($(basename(first_file_payload.name))). Second PDB ligands: $(join(second_ligands, ", ")). Check console logs."
               println("Error: " * msg)
               isfile(temp_first_path) && rm(temp_first_path, force=true); isfile(temp_second_path) && rm(temp_second_path, force=true)
               return html("<h1>Error</h1><p>$msg</p><p><a href='/'>Try again</a></p>")
           elseif isempty(second_ligands)
                msg = "No ligands found in the second PDB file ($(basename(second_file_payload.name))). First PDB ligands: $(join(first_ligands, ", ")). Check console logs."
                println("Error: " * msg)
                isfile(temp_first_path) && rm(temp_first_path, force=true); isfile(temp_second_path) && rm(temp_second_path, force=true)
                return html("<h1>Error</h1><p>$msg</p><p><a href='/'>Try again</a></p>")
           end
           
           first_ligand_options = join(["<option value='$lig'>$lig</option>" for lig in first_ligands], "")
           second_ligand_options = join(["<option value='$lig'>$lig</option>" for lig in second_ligands], "")

           # NOTE: temp_first_path and temp_second_path are passed to the /analyze route.
           # They should be cleaned up AFTER /analyze is done with them.
           # So, DO NOT delete them here if ligand selection is successful.

           html_body = """
           <!DOCTYPE html><html><head><title>juProt: Select Ligands</title>
           <link rel="icon" href="/img/favicon.ico" type="image/x-icon">
           <style>
               body { font-family: Arial, sans-serif; margin: 40px; background-color: #f4f7f6; color: #333; }
               .container { max-width: 800px; margin: auto; background-color: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
               .form-group { margin-bottom: 20px; } label { display: block; margin-bottom: 8px; font-weight: bold; }
               select { width: 100%; padding: 10px; border: 1px solid #ccc; border-radius: 4px; box-sizing: border-box;}
               button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; border-radius: 4px; font-size: 16px; }
               button:hover { background: #0056b3; } h1 { color: #0056b3; text-align: center; }
               a { color: #007bff; text-decoration: none; } a:hover { text-decoration: underline; }
               .footer-link { text-align: center; margin-top: 20px; }
           </style></head><body><div class="container">
           <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
           <h1>juProt: Select Ligands</h1>
           <form action="/analyze" method="post">
           <input type="hidden" name="first_pdb_path" value="$temp_first_path">
           <input type="hidden" name="second_pdb_path" value="$temp_second_path">
           <div class="form-group"><label for="first_ligand">First Complex Ligand (File: $(basename(first_file_payload.name))):</label>
           <select id="first_ligand" name="first_ligand" required>$first_ligand_options</select></div>
           <div class="form-group"><label for="second_ligand">Second Complex Ligand (File: $(basename(second_file_payload.name))):</label>
           <select id="second_ligand" name="second_ligand" required>$second_ligand_options</select></div>
           <button type="submit">Run Analysis</button></form>
           <div class="footer-link"><p><a href="/">Back to Home</a></p></div></div></body></html>"""
           return html(html_body)

       catch e 
           isfile(temp_first_path) && rm(temp_first_path, force=true)
           isfile(temp_second_path) && rm(temp_second_path, force=true)
           println("Error during ligand selection stage (/select-ligands route): ", sprint(showerror, e))
           Base.show_backtrace(stdout, Base.catch_backtrace())
           println()
           if e isa PythonCall.PyException
                println("--- Python Exception in /select-ligands (likely from get_ligand_names) ---")
                py_exc_obj = e.exc
                py_traceback_mod = pyimport("traceback")
                println("Python Exception (pyrepr): $(PythonCall.pyrepr(String, py_exc_obj))")
                try
                    py_traceback_mod.print_exception(pytype(py_exc_obj), py_exc_obj, py_exc_obj.__traceback__)
                catch py_print_err
                    println("Failed to print Python traceback: $py_print_err")
                end
                println("----------------------------------------------------------------------")
           end
           return html("<h1>Server Error</h1><p>An unexpected error occurred while processing PDB files for ligand selection. Details: $(sprint(showerror, e))</p><p><a href='/'>Try again</a></p>")
       end
   end

   route("/analyze", method=POST) do
       first_pdb_path = postpayload(:first_pdb_path) 
       second_pdb_path = postpayload(:second_pdb_path)
       first_ligand = postpayload(:first_ligand)
       second_ligand = postpayload(:second_ligand)

       println("Analyzing with PDBs: $first_pdb_path, $second_pdb_path and ligands: $first_ligand, $second_ligand")
       
       result = process_pdb_files(first_pdb_path, second_pdb_path, first_ligand, second_ligand)

       try
           isfile(first_pdb_path) && rm(first_pdb_path, force=true)
           isfile(second_pdb_path) && rm(second_pdb_path, force=true)
           println("Cleaned up temporary files: $first_pdb_path, $second_pdb_path")
       catch e_rm
           println("Warning: Could not clean up temporary files: $e_rm")
       end

       if haskey(result, "error")
           println("Error during analysis: ", result["error"])
           return html("<h1>Error During Analysis</h1><p>$(result["error"])</p><p><a href='/'>Try again</a></p>")
       end
       
       summary_html_safe =replace(result["summary"], "&"=>"&", "<"=>"<", ">"=>">")

       html_body = """
       <!DOCTYPE html><html><head><title>juProt Analysis Results</title>
       <style>
           body { font-family: Arial, sans-serif; margin: 40px; background-color: #f4f7f6; color: #333; }
           .container { max-width: 900px; margin: auto; background-color: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
           .result-section { margin-top: 30px; padding-bottom: 20px; border-bottom: 1px solid #eee; }
           .result-section:last-child { border-bottom: none; }
           h1 { color: #0056b3; text-align: center; margin-bottom: 30px;} h2 { color: #007bff; margin-bottom: 15px; }
           a { color: #007bff; text-decoration: none; } a:hover { text-decoration: underline; }
           img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; margin-top: 10px; }
           pre { background: #e9ecef; padding: 15px; border-radius: 5px; white-space: pre-wrap; word-wrap: break-word; font-size: 0.9em; }
           ul { list-style-type: none; padding-left: 0; } ul li { margin-bottom: 8px; }
           ul li a { display: inline-block; padding: 8px 12px; background-color: #007bff; color: white; border-radius: 4px; }
           ul li a:hover { background-color: #0056b3; } .footer-link { text-align: center; margin-top: 30px; }
       </style></head><body><div class="container">
       <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
       <h1>Analysis Results</h1>
       <p>Comparison of H-bond interactions for the selected ligands in the First and Second Complexes.</p>
       <div class="result-section"><h2>Downloads</h2><ul>
       <li><a href="$(result["comparison_table"])" download>Download Comparison Table (CSV)</a></li>
       <li><a href="$(result["detailed_interactions"])" download>Download Detailed Interactions (CSV)</a></li>
       <li><a href="$(result["plot_file"])" download>Download Bar Chart (PNG)</a></li>
       </ul></div><div class="result-section"><h2>Residue Interaction Plot</h2>
       <img src="$(result["plot_file"])" alt="Residue Interactions Bar Chart"></div>
       <div class="result-section"><h2>Analytical Summary</h2><pre>$(summary_html_safe)</pre></div>
       <div class="footer-link"><p><a href="/">Run Another Analysis</a></p></div></div></body></html>"""
       return html(html_body)
   end

   route("/how-to-use") do
       html("""
       <!DOCTYPE html>
<html>
<head>
    <title>juProt: How to Use & Applications</title>
    <link rel="icon" href="/img/favicon.ico" type="image/x-icon">
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background-color: #f4f7f6; color: #333; }
        .container { max-width: 800px; margin: auto; background-color: #fff; padding: 20px 40px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
        h1, h2, h3 { color: #0056b3; }
        h1 { text-align: center; margin-bottom: 30px; }
        h2 { margin-top: 25px; border-bottom: 1px solid #eee; padding-bottom: 5px;}
        h3 { margin-top: 20px; color: #007bff; }
        p, li { line-height: 1.6; }
        a { color: #007bff; text-decoration: none; }
        a:hover { text-decoration: underline; }
        ol, ul { padding-left: 20px; }
        .footer-link { text-align: center; margin-top: 30px; }
        code { background-color: #e9ecef; padding: 2px 4px; border-radius: 3px; font-family: monospace; }
    </style>
</head>
<body>
    <div class="container">
    <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
        <h1>How to Use juProt & Its Applications</h1>

        <h2>Protocol</h2>
        <ol>
            <li><strong>Prepare PDB Files</strong>: Obtain two protein-ligand complex PDB files (e.g., from RCSB PDB or molecular modeling). Ensure HETATM records for ligands are present.</li>
            <li><strong>Access the App</strong>: Visit juProt at <code>https://juprot.info/</code>.</li>
            <li><strong>Upload Files</strong>:
                <ul>
                    <li>Upload the first PDB file under "First Complex PDB File".</li>
                    <li>Upload the second PDB file under "Second Complex PDB File".</li>
                    <li>Click "Load Ligands".</li>
                </ul>
            </li>
            <li><strong>Select Ligands</strong>:
                <ul>
                    <li>juProt will auto-detect potential ligands from each PDB.</li>
                    <li>Select the specific ligand of interest from the dropdown menu for each complex.</li>
                    <li>Click "Run Analysis".</li>
                </ul>
            </li>
            <li><strong>View Results</strong>:
                <ul>
                    <li>Examine the textual "Analytical Summary" for a quick overview.</li>
                    <li>View the "Residue Interaction Plot" (bar chart) for a visual comparison of H-bond frequencies per residue.</li>
                    <li>Download the "Comparison Table (CSV)" for a structured summary of differences and commonalities.</li>
                    <li>Download "Detailed Interactions (CSV)" for a complete list of all H-bonds for both complexes with their geometric parameters.</li>
                </ul>
            </li>
        </ol>

        <h2>Applications & Use Cases for juProt</h2>
        <p>juProt is designed to provide rapid comparative insights into protein-ligand interactions, with a current focus on hydrogen bonds. Here are some scenarios where juProt can be particularly useful:</p>

        <h3>1. Understanding the Impact of Mutations</h3>
        <p><strong>Scenario:</strong> You have a wild-type protein-ligand structure and a mutant form (e.g., from a SNP or site-directed mutagenesis) bound to the same ligand.</p>
        <p><strong>How juProt Helps:</strong> Compare the native-ligand and mutant-ligand complexes. juProt highlights how the mutation alters the H-bond network, which can help explain changes in binding affinity, drug efficacy, or resistance mechanisms. (e.g., comparing a wild-type kinase-inhibitor complex with a gatekeeper mutant-inhibitor complex).</p>

        <h3>2. Comparing Different Ligands to the Same Target</h3>
        <p><strong>Scenario:</strong> You have several drug candidates or chemical probes binding to the same protein target.</p>
        <p><strong>How juProt Helps:</strong> Compare Protein+LigandA with Protein+LigandB. juProt helps identify which ligand forms more/different H-bonds and which residues are key common or unique H-bond partners, aiding in SAR studies and lead optimization.</p>

        <h3>3. Analyzing Ligand Binding to Different Protein Conformations or Isoforms</h3>
        <p><strong>Scenario:</strong> A protein exists in different states (e.g., active/inactive) or as different isoforms, and you have structures of a ligand bound to these variants.</p>
        <p><strong>How juProt Helps:</strong> Compare ProteinState1+Ligand with ProteinState2+Ligand. juProt can reveal how protein structural changes influence H-bonding with a common ligand.</p>

        <h3>4. Validating Docking Poses</h3>
        <p><strong>Scenario:</strong> You have multiple potential binding poses for a ligand from molecular docking.</p>
        <p><strong>How juProt Helps:</strong> Compare the H-bond profile of different docked poses or a docked pose against an experimental structure (if available) to assess H-bonding consistency.</p>

        <h3>5. Educational Purposes</h3>
        <p><strong>Scenario:</strong> Teaching students about protein-ligand interactions.</p>
        <p><strong>How juProt Helps:</strong> Provides an easy-to-use tool for students to explore H-bonds and understand differential interactions without complex software or scripting.</p>

        <h3>Future Enhancements</h3>
        <p>juProt is an actively developed open-source project. Future versions aim to expand comparative analysis to other interaction types (hydrophobic contacts, π-stacking, salt bridges, etc.) for a more comprehensive view of protein-ligand interactomes.</p>

        <div class="footer-link"><p><a href="/">Back to Home</a></p></div>
    </div>
</body>
</html>""")
   end

   route("/about") do
    html("""
    <!DOCTYPE html>
    <html>
    <head>
        <title>About juProt</title>
        <link rel="icon" href="/img/favicon.ico" type="image/x-icon">
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; background-color: #f4f7f6; color: #333; }
            .container { max-width: 800px; margin: auto; background-color: #fff; padding: 20px 40px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
            h1, h2, h3 { color: #0056b3; }
            h1 { text-align: center; margin-bottom: 30px; }
            h2 { margin-top: 25px; border-bottom: 1px solid #eee; padding-bottom: 5px;}
            h3 { margin-top: 20px; color: #007bff; }
            p, li { line-height: 1.6; }
            a { color: #007bff; text-decoration: none; }
            a:hover { text-decoration: underline; }
            .footer-link { text-align: center; margin-top: 30px; }
            code { background-color: #e9ecef; padding: 2px 4px; border-radius: 3px; font-family: monospace; }
        </style>
    </head>
    <body>
        <div class="container">
        <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
            <h1>About juProt</h1>

            <p>juProt (<code>https://juprot.info/</code>) is an open-source web application designed to facilitate the rapid and user-friendly comparative analysis of protein-ligand interaction networks, with an initial focus on hydrogen bonds. Understanding how ligands interact with their protein targets, and how these interactions change due to mutations or when comparing different ligands, is fundamental in structural biology, bioinformatics, and drug discovery.</p>

            <h2>Motivation</h2>
            <p>While powerful tools exist for analyzing interactions in a single protein-ligand complex, comparing these interactions across two different complexes often requires manual data extraction, scripting, and collation of results. juProt aims to simplify this process, making comparative interaction analysis accessible to a broader range of researchers and students without requiring extensive computational expertise.</p>

            <h2>Core Technology</h2>
            
            <p>juProt is developed using the Julia programming language, leveraging the high-performance Genie.jl web framework for its backend and user interface. The core interaction detection is powered by the well-established Protein-Ligand Interaction Profiler (PLIP), a Python-based tool. juProt interfaces with PLIP using the PythonCall.jl package, allowing seamless integration of PLIP's robust algorithms.</p>

            <h2>Current Features</h2>
            <ul>
                <li>Upload of two PDB files for comparison.</li>
                <li>Automated identification of potential ligands with user selection.</li>
                <li>Detection and quantification of hydrogen bonds for each complex.</li>
                <li>Generation of:
                    <ul>
                        <li>A comparative summary table (CSV) highlighting common and differential H-bonds and interacting residues.</li>
                        <li>A detailed list of all H-bonds for both complexes (CSV).</li>
                        <li>A bar chart visually comparing H-bond counts per residue across the two complexes (PNG).</li>
                        <li>An on-page analytical summary of key findings.</li>
                    </ul>
                </li>
            </ul>

            <h2>Future Development</h2>
            <p>juProt is an ongoing project. Future development plans include:</p>
            <ul>
                <li>Expanding comparative analysis to other interaction types (hydrophobic, pi-stacking, salt bridges, etc.).</li>
                <li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li>
                <li>Allowing analysis of more than two complexes.</li>
                <li>User-configurable parameters for interaction detection.</li>
            </ul>

            <h2>Development Team</h2>
            <p>juProt was conceived and developed by:<br>
               <p><a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><p>
               <p>Assistant Professor, Department of Biotechnology, Sri Ramachandra Institute of Higher Education and Research, Chennai - 600116, Tamil Nadu, India.</p>
                
                <p>Deepak S P, MSc Biotechnology</p>
                <p>Siva V, MSc Biotechnology</p>
                <p>Surya Sekaran, [PhD]</p>
            
            <p>We also acknowledge the developers of the core libraries used in juProt, including Julia, Genie.jl, PythonCall.jl, PLIP, OpenBabel, and Plots.jl.</p>

            <h2>Open Source & Citation</h2>
            <p>juProt is an open-source project. The source code is available on GitHub at <a href="https://github.com/drbenedictpaul/juprot" target="_blank">https://github.com/drbenedictpaul/juprot</a>.</p>
            <p>We encourage contributions and feedback from the community.</p>
            <p>If you use juProt in your research, please cite:<br>
               <em>[Manuscript is in communication. For now, please cite our GitHub repository.]</em></p>

            <h2>Contact/Feedback</h2>
            <p>For questions, suggestions, or to report issues, please visit our GitHub issues page at <a href="https://github.com/drbenedictpaul/juprot/issues" target="_blank">GitHub Issues</a> or contact <a href="mailto:benedictpaulc@sriramachandra.edu.in">benedictpaulc@sriramachandra.edu.in</a>.</p>

            <div class="footer-link"><p><a href="/">Back to Home</a></p></div>
        </div>
    </body>
    </html>
    """)
end

end # module JuProtGUI