module JuProtGUI
   using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
   using CSV, DataFrames, Plots, PythonCall

   include("../lib/output_results.jl")
   include("../lib/detect_hbonds.jl")

   function get_ligand_names(pdb_file)
       try
           plip = pyimport("plip.structure.preparation")
           mol = plip.PDBComplex()
           mol.load_pdb(pdb_file)
           ligands = Set{String}()
           for res in mol.ligands
               push!(ligands, res.resn)
           end
           println("Detected potential ligands: ", join(ligands, ", "))
           return collect(ligands)
       catch e
           println("Error reading PDB file ", pdb_file, ": ", sprint(showerror, e))
           return []
       end
   end

   function process_pdb_files(first_pdb_path, second_pdb_path, first_ligand, second_ligand)
       println("Starting process_pdb_files with files: ", first_pdb_path, ", ", second_pdb_path)
       complex_results = Dict{String, Dict}()
       output_dir = "public/outputs"
       mkpath(output_dir)
       first_path = joinpath(output_dir, "first_complex.pdb")
       second_path = joinpath(output_dir, "second_complex.pdb")
       cp(first_pdb_path, first_path, force=true)
       cp(second_pdb_path, second_path, force=true)
       pdb_files = [first_path, second_path]
       ligand_resnames = [first_ligand, second_ligand]

       for (i, pdb_file) in enumerate(pdb_files)
           try
               println("Processing PDB file: ", pdb_file)
               ligand_resname = ligand_resnames[i]
               hbonds = detect_hbonds(pdb_file, ligand_resname)
               close_contacts = hbonds
               println("Found ", length(close_contacts), " H-bonds")
               interacting_residues = Dict{String, Int}()
               for (p_atom, _, _, _) in close_contacts
                   res_key = "$(p_atom.resn) $(p_atom.resi)"
                   interacting_residues[res_key] = get(interacting_residues, res_key, 0) + 1
               end
               max_residue = ""
               max_count = 0
               for (res, count) in interacting_residues
                   if count > max_count
                       max_residue = res
                       max_count = count
                   end
               end
               complex_name = basename(pdb_file)
               complex_results[complex_name] = Dict(
                   :close_contacts => close_contacts,
                   :hbonds => hbonds,
                   :interacting_residues => interacting_residues,
                   :max_residue => max_residue,
                   :max_count => max_count,
                   :ligand_resname => ligand_resname
               )
           catch e
               println("Error processing PDB file ", pdb_file, ": ", sprint(showerror, e))
               return Dict("error" => "Failed to process $pdb_file: $(sprint(showerror, e))")
           end
       end
       output_file = joinpath(output_dir, "comparison_table.csv")
       plot_file = joinpath(output_dir, "residue_interactions.png")
       try
           println("Generating comparison and summary")
           compare_and_save_results(complex_results, output_file, plot_file)
           summary = capture_analytical_summary(complex_results)
           println("Captured summary: ", summary)
           return Dict(
               "comparison_table" => output_file,
               "detailed_interactions" => joinpath(output_dir, "detailed_interactions.csv"),
               "plot_file" => plot_file,
               "summary" => summary
           )
       catch e
           println("Error during analysis: ", sprint(showerror, e))
           return Dict("error" => "Analysis failed: $(sprint(showerror, e))")
       end
   end

   function capture_analytical_summary(complex_results)
       println("Capturing analytical summary")
       try
           output = Pipe()
           redirect_stdout(output) do
               print_analytical_summary(complex_results)
           end
           close(output.in)
           summary = String(read(output))
           if isempty(summary)
               println("Warning: Empty analytical summary captured")
           else
               println("Summary captured successfully: ", summary)
           end
           return summary
       catch e
           println("Error capturing summary: ", sprint(showerror, e))
           return "Error capturing summary: $(sprint(showerror, e))"
       end
   end

   route("/") do
       html("""
       <!DOCTYPE html>
       <html>
       <head>
           <title>juProt</title>
           <style>
               body { font-family: Arial, sans-serif; margin: 40px; }
               .container { max-width: 800px; margin: auto; }
               .form-group { margin-bottom: 20px; }
               label { display: block; margin-bottom: 5px; }
               input[type=file] { width: 100%; padding: 8px; }
               button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; }
               button:hover { background: #0056b3; }
               a { color: #007bff; text-decoration: none; }
               a:hover { text-decoration: underline; }
           </style>
       </head>
       <body>
           <div class="container">
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
               <p><a href="/how-to-use">How to Use</a></p>
           </div>
       </body>
       </html>
       """)
   end

   route("/select-ligands", method=POST) do
       first_file = filespayload("first_pdb")
       second_file = filespayload("second_pdb")
       output_dir = "public/outputs"
       mkpath(output_dir)
       first_path = joinpath(output_dir, "temp_first.pdb")
       second_path = joinpath(output_dir, "temp_second.pdb")
       write(first_path, first_file.data)
       write(second_path, second_file.data)

       first_ligands = get_ligand_names(first_path)
       second_ligands = get_ligand_names(second_path)

       if isempty(first_ligands) || isempty(second_ligands)
           println("Error: No ligands found in one or both PDB files: ", first_ligands, ", ", second_ligands)
           return html("<h1>Error</h1><p>No ligands found in one or both PDB files.</p>")
       end

       html("""
       <!DOCTYPE html>
       <html>
       <head>
           <title>juProt: Select Ligands</title>
           <style>
               body { font-family: Arial, sans-serif; margin: 40px; }
               .container { max-width: 800px; margin: auto; }
               .form-group { margin-bottom: 20px; }
               label { display: block; margin-bottom: 5px; }
               select { width: 100%; padding: 8px; }
               button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; }
               button:hover { background: #0056b3; }
               a { color: #007bff; text-decoration: none; }
               a:hover { text-decoration: underline; }
           </style>
       </head>
       <body>
           <div class="container">
               <h1>juProt: Select Ligands</h1>
               <form action="/analyze" method="post">
                   <input type="hidden" name="first_pdb_path" value="$first_path">
                   <input type="hidden" name="second_pdb_path" value="$second_path">
                   <div class="form-group">
                       <label for="first_ligand">First Complex Ligand:</label>
                       <select id="first_ligand" name="first_ligand" required>
                           $(join(["<option value='$ligand'>$ligand</option>" for ligand in first_ligands], ""))
                       </select>
                   </div>
                   <div class="form-group">
                       <label for="second_ligand">Second Complex Ligand:</label>
                       <select id="second_ligand" name="second_ligand" required>
                           $(join(["<option value='$ligand'>$ligand</option>" for ligand in second_ligands], ""))
                       </select>
                   </div>
                   <button type="submit">Run Analysis</button>
               </form>
               <p><a href="/">Back to Home</a></p>
           </div>
       </body>
       </html>
       """)
   end

   route("/analyze", method=POST) do
       first_pdb_path = postpayload(:first_pdb_path)
       second_pdb_path = postpayload(:second_pdb_path)
       first_ligand = postpayload(:first_ligand)
       second_ligand = postpayload(:second_ligand)
       println("Analyzing with ligands: ", first_ligand, ", ", second_ligand)
       result = process_pdb_files(first_pdb_path, second_pdb_path, first_ligand, second_ligand)
       if haskey(result, "error")
           println("Error during analysis: ", result["error"])
           return html("<h1>Error</h1><p>$(result["error"])</p>")
       end
       html("""
       <!DOCTYPE html>
       <html>
       <head>
           <title>juProt Analysis Results</title>
           <style>
               body { font-family: Arial, sans-serif; margin: 40px; }
               .container { max-width: 800px; margin: auto; }
               .result { margin-top: 20px; }
               a { color: #007bff; text-decoration: none; }
               a:hover { text-decoration: underline; }
               img { max-width: 100%; height: auto; }
               pre { background: #f8f9fa; padding: 10px; border-radius: 5px; }
           </style>
       </head>
       <body>
           <div class="container">
               <h1>Analysis Results</h1>
               <p>Comparison of H-bond interactions for the First and Second Complexes:</p>
               <div class="result">
                   <h2>Downloads</h2>
                   <p><a href="/outputs/comparison_table.csv">Download Comparison Table (CSV)</a></p>
                   <p><a href="/outputs/detailed_interactions.csv">Download Detailed Interactions (CSV)</a></p>
                   <p><a href="/outputs/residue_interactions.png">Download Bar Chart (PNG)</a></p>
               </div>
               <div class="result">
                   <h2>Residue Interaction Plot</h2>
                   <img src="/outputs/residue_interactions.png" alt="Residue Interactions">
               </div>
               <div class="result">
                   <h2>Analytical Summary</h2>
                   <pre>$(result["summary"])</pre>
               </div>
               <p><a href="/">Run Another Analysis</a></p>
           </div>
       </body>
       </html>
       """)
   end

   route("/how-to-use") do
       html("""
       <!DOCTYPE html>
       <html>
       <head>
           <title>juProt: How to Use</title>
           <style>
               body { font-family: Arial, sans-serif; margin: 40px; }
               .container { max-width: 800px; margin: auto; }
               h1, h2 { color: #007bff; }
               p { line-height: 1.6; }
               a { color: #007bff; text-decoration: none; }
               a:hover { text-decoration: underline; }
           </style>
       </head>
       <body>
           <div class="container">
               <h1>How to Use juProt</h1>
               <p>juProt is a web-based tool for comparing protein-ligand hydrogen bond interactions between two complexes.</p>
               <h2>Protocol</h2>
               <ol>
                   <li><strong>Prepare PDB Files</strong>: Obtain two protein-ligand complex PDB files (e.g., from PDB).</li>
                   <li><strong>Access the App</strong>: Visit <a href="/">juProt</a>.</li>
                   <li><strong>Upload Files</strong>:
                       <ul>
                           <li>Upload the first PDB file under "First Complex PDB File".</li>
                           <li>Upload the second PDB file under "Second Complex PDB File".</li>
                           <li>Select ligands from the dropdown menus for each complex.</li>
                       </ul>
                   </li>
                   <li><strong>Run Analysis</strong>: Click "Run Analysis" to compare hydrogen bonds.</li>
                   <li><strong>View Results</strong>:
                       <ul>
                           <li>Download the comparison table (CSV) for interaction counts.</li>
                           <li>Download detailed interactions (CSV) for specific interactions.</li>
                           <li>View or download the bar chart (PNG) comparing residue interactions.</li>
                           <li>Read the analytical summary for insights.</li>
                       </ul>
                   </li>
               </ol>
               <h2>Applications</h2>
               <ul>
                   <li><strong>Structural Biology</strong>: Compare protein-ligand complexes.</li>
                   <li><strong>Drug Design</strong>: Analyze ligand binding differences.</li>
                   <li><strong>Protein Engineering</strong>: Assess interaction changes.</li>
               </ul>
               <p><a href="/">Back to Home</a></p>
           </div>
       </body>
       </html>
       """)
   end

   end