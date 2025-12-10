module JuProtGUI
   using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
   using CSV, DataFrames, PythonCall, Printf, Downloads
   
   const LIB_PATH = joinpath(@__DIR__, "..", "lib")
   const OUTPUT_DIR = joinpath(@__DIR__, "..", "public", "outputs")
   
   include(joinpath(LIB_PATH, "output_results.jl"))
   include(joinpath(LIB_PATH, "detect_hbonds.jl"))

   # --- PyMOL Script Generator ---
   function generate_pymol_script(complex_results)
       script_content = IOBuffer()
       path1 = abspath(complex_results["first_complex.pdb"][:original_path]); path2 = abspath(complex_results["second_complex.pdb"][:original_path])
       name1 = "complex_1_" * replace(basename(path1), r"[^a-zA-Z0-9]" => "_"); name2 = "complex_2_" * replace(basename(path2), r"[^a-zA-Z0-9]" => "_")
       println(script_content, "load \"$path1\", $name1"); println(script_content, "load \"$path2\", $name2"); println(script_content, "align $name2, $name1")
       println(script_content, "bg_color white"); println(script_content, "hide everything, all"); println(script_content, "show cartoon, polymer.protein"); println(script_content, "color gray80, polymer.protein")
       function get_resi_string(res_dict)
           selections = String[]; for res_key in keys(res_dict); parts = split(res_key); if length(parts) >= 2; push!(selections, "resi $(parts[2])"); end; end; return join(selections, " or ")
       end
       res1_str = get_resi_string(complex_results["first_complex.pdb"][:interacting_residues]); res2_str = get_resi_string(complex_results["second_complex.pdb"][:interacting_residues])
       if !isempty(res1_str); println(script_content, "select resi_c1, ($res1_str) and $name1"); println(script_content, "color marine, resi_c1"); println(script_content, "show sticks, resi_c1"); println(script_content, "util.cbac('blue', 'resi_c1')"); end
       if !isempty(res2_str); println(script_content, "select resi_c2, ($res2_str) and $name2"); println(script_content, "color firebrick, resi_c2"); println(script_content, "show sticks, resi_c2"); println(script_content, "util.cba(214, 'resi_c2')"); end
       lig1 = complex_results["first_complex.pdb"][:ligand_resname]; lig2 = complex_results["second_complex.pdb"][:ligand_resname]
       println(script_content, "show sticks, (resn $lig1 and $name1) or (resn $lig2 and $name2)"); println(script_content, "zoom (resn $lig1 and $name1) or (resn $lig2 and $name2)")
       return String(take!(script_content))
   end

   # --- Backend Logic ---
   function get_ligand_names(pdb_file_path::String); try; plip = pyimport("plip.structure.preparation"); mol = plip.PDBComplex(); mol.load_pdb(pdb_file_path); ligands = Set{String}(); if PythonCall.pyhasattr(mol, "ligands") && mol.ligands !== PythonCall.pybuiltins.None; for res in mol.ligands; if PythonCall.pyhasattr(res, "hetid"); name_py = res.hetid; if name_py !== PythonCall.pybuiltins.None; name = pyconvert(String, name_py); if !isempty(name) && name != "HOH"; push!(ligands, name); end; end; end; end; end; if isempty(ligands) && PythonCall.pyhasattr(mol, "res_het"); for res in mol.res_het; if PythonCall.pyhasattr(res, "hetid"); name = pyconvert(String, res.hetid); if !isempty(name) && name != "HOH"; push!(ligands, name); end; end; end; end; return collect(ligands); catch e; return String[]; end; end
   function process_pdb_files(first_pdb_path::String, second_pdb_path::String, first_ligand::String, second_ligand::String); complex_results = Dict{String, Dict}(); mkpath(OUTPUT_DIR); first_complex_key = "first_complex.pdb"; second_complex_key = "second_complex.pdb"; pdb_files_to_process = [(path=first_pdb_path, ligand=first_ligand, key=first_complex_key), (path=second_pdb_path, ligand=second_ligand, key=second_complex_key)]; all_processed_successfully = true; for item in pdb_files_to_process; try; interactions = detect_hbonds(item.path, item.ligand); interacting_residues = Dict{String, Int}(); for detail in interactions; res_key = "$(detail.prot_resn) $(detail.prot_resi)"; interacting_residues[res_key] = get(interacting_residues, res_key, 0) + 1; end; complex_results[item.key] = Dict(:hbonds => interactions, :interacting_residues => interacting_residues, :ligand_resname => item.ligand, :original_path => item.path); catch e; all_processed_successfully = false; end; end; if !all_processed_successfully; return Dict("error" => "Processing failed."); end; output_comparison_csv = joinpath(OUTPUT_DIR, "comparison_table.csv"); plot_file_path = joinpath(OUTPUT_DIR, "residue_interactions.png"); detailed_interactions_csv = joinpath(OUTPUT_DIR, "detailed_interactions.csv"); try; compare_and_save_results(complex_results, first_complex_key, second_complex_key, output_comparison_csv, detailed_interactions_csv, plot_file_path); summary_text = capture_analytical_summary(complex_results, first_complex_key, second_complex_key); return Dict("comparison_table" => "/outputs/comparison_table.csv", "detailed_interactions" => "/outputs/detailed_interactions.csv", "plot_file" => "/outputs/residue_interactions.png", "summary" => summary_text, "session_data" => complex_results); catch e; return Dict("error" => "Visualization failed: $(sprint(showerror, e))"); end; end
   function capture_analytical_summary(complex_results, first_key, second_key); try; io=IOBuffer(); print_analytical_summary(io,complex_results,first_key,second_key); return String(take!(io)); catch; return "Summary failed."; end; end

   # --- Routes ---

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
            .form-group { margin-bottom: 20px; position: relative; }
            label { display: block; margin-bottom: 8px; font-weight: bold; }
            input[type=file], input[type=text] { width: calc(100% - 22px); padding: 10px; border: 1px solid #ccc; border-radius: 4px; }
            button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; border-radius: 4px; font-size: 16px; }
            button:hover { background: #0056b3; }
            h1 { color: #0056b3; text-align: center; }
            a { color: #007bff; text-decoration: none; }
            a:hover { text-decoration: underline; }
            .footer-links-container { display: flex; justify-content: center; gap: 30px; margin-top: 40px; padding-top: 20px; border-top: 1px solid #eee; }
            .footer-links-container p { margin: 0; }
            .divider { text-align: center; font-weight: bold; color: #aaa; margin: 20px 0; }
            .validation-popup {
                display: none; 
                position: absolute;
                top: 35px;
                left: 110px;
                background-color: #fff;
                border: 1px solid #ccc;
                border-radius: 4px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.15);
                padding: 8px 12px;
                z-index: 100;
                white-space: nowrap;
                font-size: 14px;
                color: #333;
                align-items: center;
            }
            .validation-popup::before { content: ''; position: absolute; bottom: 100%; left: 20px; border-width: 7px; border-style: solid; border-color: transparent transparent #ccc transparent; }
            .validation-popup::after { content: ''; position: absolute; bottom: 100%; left: 21px; border-width: 6px; border-style: solid; border-color: transparent transparent #fff transparent; }
            .validation-popup .icon {
                display: inline-block;
                background-color: #ff8552;
                color: white;
                width: 16px;
                height: 16px;
                border-radius: 3px;
                text-align: center;
                font-weight: bold;
                line-height: 16px;
                margin-right: 8px;
                font-size: 12px;
            }
        </style>
        <script>
            function validateForm() {
                document.getElementById('alert_1').style.display = 'none';
                document.getElementById('alert_2').style.display = 'none';
                var file1 = document.getElementById('first_pdb').value;
                var id1 = document.getElementById('first_pdb_id').value.trim();
                var complex1_provided = (file1 !== "" || id1 !== "");
                var file2 = document.getElementById('second_pdb').value;
                var id2 = document.getElementById('second_pdb_id').value.trim();
                var complex2_provided = (file2 !== "" || id2 !== "");
                var isValid = true;
                if (!complex1_provided) {
                    document.getElementById('alert_1').style.display = 'flex';
                    isValid = false;
                }
                if (!complex2_provided) {
                    document.getElementById('alert_2').style.display = 'flex';
                    isValid = false;
                }
                return isValid;
            }
        </script>
    </head>
    <body>
        <div class="container">
            <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
            <h1>juProt: Protein-Ligand Interaction Analyzer</h1>
            <!-- <p style="text-align: center; color: #666;">Comparative analysis of the complete protein-ligand interactome.</p> -->
            <form action="/select-ligands" method="post" enctype="multipart/form-data" onsubmit="return validateForm();">
                <div class="form-group">
                    <label for="first_pdb">Upload First Complex (PDB File):</label>
                    <input type="file" id="first_pdb" name="first_pdb" accept=".pdb">
                    <div id="alert_1" class="validation-popup"><span class="icon">!</span>Please select a file or enter a PDB ID.</div>
                </div>
                <div class="form-group">
                    <label for="second_pdb">Upload Second Complex (PDB File):</label>
                    <input type="file" id="second_pdb" name="second_pdb" accept=".pdb">
                    <div id="alert_2" class="validation-popup"><span class="icon">!</span>Please select a file or enter a PDB ID.</div>
                </div>
                <div class="divider">OR</div>
                <div class="form-group">
                    <label for="first_pdb_id">Enter First PDB ID (e.g., 3EQM):</label>
                    <input type="text" id="first_pdb_id" name="first_pdb_id" placeholder="PDB ID">
                </div>
                <div class="form-group">
                    <label for="second_pdb_id">Enter Second PDB ID:</label>
                    <input type="text" id="second_pdb_id" name="second_pdb_id" placeholder="PDB ID">
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
       first_file_payload = haskey(filespayload(), :first_pdb) ? filespayload(:first_pdb) : nothing
       second_file_payload = haskey(filespayload(), :second_pdb) ? filespayload(:second_pdb) : nothing
       first_pdb_id = postpayload(:first_pdb_id, ""); second_pdb_id = postpayload(:second_pdb_id, "")
       mkpath(OUTPUT_DIR); temp_first_path = ""; temp_second_path = ""; first_filename = ""; second_filename = ""
       try
           if !isempty(first_pdb_id); temp_first_path = joinpath(OUTPUT_DIR, "$(first_pdb_id).pdb"); Downloads.download("https://files.rcsb.org/download/$(first_pdb_id).pdb", temp_first_path); first_filename = "$(first_pdb_id).pdb"; elseif first_file_payload !== nothing && !isempty(first_file_payload.data); temp_first_path = joinpath(OUTPUT_DIR, "temp_first_" * string(randn())[3:8] * ".pdb"); write(temp_first_path, first_file_payload.data); first_filename = first_file_payload.name; else; return html("<h1>Error</h1><p>No input for the first complex.</p>"); end
           if !isempty(second_pdb_id); temp_second_path = joinpath(OUTPUT_DIR, "$(second_pdb_id).pdb"); Downloads.download("https://files.rcsb.org/download/$(second_pdb_id).pdb", temp_second_path); second_filename = "$(second_pdb_id).pdb"; elseif second_file_payload !== nothing && !isempty(second_file_payload.data); temp_second_path = joinpath(OUTPUT_DIR, "temp_second_" * string(randn())[3:8] * ".pdb"); write(temp_second_path, second_file_payload.data); second_filename = second_file_payload.name; else; return html("<h1>Error</h1><p>No input for the second complex.</p>"); end
           first_ligands = get_ligand_names(temp_first_path); second_ligands = get_ligand_names(temp_second_path)
           if isempty(first_ligands) && isempty(second_ligands); return html("<h1>Error</h1><p>No ligands found.</p>"); end
           first_ligand_options = join(["<option value='$l'>$l</option>" for l in first_ligands], ""); second_ligand_options = join(["<option value='$l'>$l</option>" for l in second_ligands], "")
           html_body = """<!DOCTYPE html><html><head><title>juProt: Select Ligands</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.form-group{margin-bottom:20px}label{display:block;margin-bottom:8px;font-weight:700}select{width:100%;padding:10px;border:1px solid #ccc;border-radius:4px;box-sizing:border-box}button{padding:10px 20px;background:#007bff;color:#fff;border:none;cursor:pointer;border-radius:4px;font-size:16px}button:hover{background:#0056b3}h1{color:#0056b3;text-align:center}.footer-link{text-align:center;margin-top:20px}a{color:#007bff;text-decoration:none}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>Select Target Ligands</h1><form action=/analyze method=post><input type=hidden name=first_pdb_path value=$temp_first_path><input type=hidden name=second_pdb_path value=$temp_second_path><div class=form-group><label>First Complex Ligand ($first_filename):</label><select name=first_ligand required>$first_ligand_options</select></div><div class=form-group><label>Second Complex Ligand ($second_filename):</label><select name=second_ligand required>$second_ligand_options</select></div><button type=submit>Run Analysis</button></form><div class=footer-link><p><a href=/>Back</a></p></div></div></body></html>"""
           return html(html_body)
       catch e; return html("<h1>Error</h1><p>Failed to process PDBs. Details: $(sprint(showerror,e))</p>"); end
   end

   SESSION_DATA = Ref{Any}(nothing)

   route("/analyze", method=POST) do
       first_pdb_path=postpayload(:first_pdb_path); second_pdb_path=postpayload(:second_pdb_path); first_ligand=postpayload(:first_ligand); second_ligand=postpayload(:second_ligand)
       result = process_pdb_files(first_pdb_path, second_pdb_path, first_ligand, second_ligand)
       if haskey(result, "session_data"); SESSION_DATA[] = result["session_data"]; end
       if haskey(result,"error"); return html("<h1>Error</h1><p>$(result["error"])</p>"); end
       summary_html_safe = replace(result["summary"], "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;", "\n"=>"<br>")
       html_body="""<!DOCTYPE html><html><head><title>juProt Analysis Results</title><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:950px;margin:auto;background-color:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.result-section{margin-top:20px;padding-bottom:20px;border-bottom:1px solid #eee}.result-section:last-child{border-bottom:none}h1{color:#0056b3;text-align:center;margin-bottom:30px}h2{color:#007bff;margin-bottom:5px;margin-top:0}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}img{max-width:100%;height:auto;border:none;margin-top:0;display:block;margin-left:auto;margin-right:auto}pre{background:#e9ecef;padding:15px;border-radius:5px;white-space:pre-wrap;word-wrap:break-word;font-size:.9em}ul{list-style-type:none;padding-left:0}ul li{margin-bottom:8px}ul li a{display:inline-block;padding:8px 12px;background-color:#007bff;color:#fff;border-radius:4px}ul li a:hover{background-color:#0056b3}.footer-link{text-align:center;margin-top:30px}</style></head><body><div class=container><div style=text-align:center;margin-bottom:20px><img src=/img/juProt_logo.png alt="JuProt Logo" style=max-height:70px;margin-top:10px></div><h1>Analysis Results</h1><p>Comparison of full interaction profiles for the selected ligands.</p><div class=result-section><h2>Downloads</h2><ul><li><a href="$(result["comparison_table"])" download>Download Summary CSV</a></li><li><a href="$(result["detailed_interactions"])" download>Download Detailed CSV</a></li><li><a href="$(result["plot_file"])" download>Download Chart Image (PNG)</a></li><li><a href="/pymol-script" download="juprot_session.pml">Download PyMOL Session (.pml)</a></li></ul></div><div class=result-section><h2>Residue Interaction Profile</h2><img src="$(result["plot_file"])" alt="Interaction Chart"></div><div class=result-section><h2>Analytical Summary</h2><pre>$(summary_html_safe)</pre></div><div class=footer-link><p><a href=/>Run Another Analysis</a></p></div></div></body></html>"""
       return html(html_body)
   end

   route("/pymol-script") do; if SESSION_DATA[]!==nothing; script_text=generate_pymol_script(SESSION_DATA[]); SESSION_DATA[]=nothing; return script_text; else; return"No session data."; end; end

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
                    <li>View the "Residue Interaction Plot" (bar chart) for a visual comparison of interaction frequencies per residue.</li>
                    <li>Download the "Comparison Table (CSV)" for a structured summary of differences and commonalities.</li>
                    <li>Download "Detailed Interactions (CSV)" for a complete list of all interactions for both complexes with their geometric parameters.</li>
                </ul>
            </li>
        </ol>

        <h2>Applications & Use Cases for juProt</h2>
        <p>juProt is designed to provide rapid comparative insights into protein-ligand interactions, with a comprehensive focus on the full interactome (Hydrogen bonds, Hydrophobic contacts, Salt Bridges, Pi-Stacking, etc.).</p>

        <h3>1. Understanding the Impact of Mutations</h3>
        <p><strong>Scenario:</strong> You have a wild-type protein-ligand structure and a mutant form (e.g., from a SNP or site-directed mutagenesis) bound to the same ligand.</p>
        <p><strong>How juProt Helps:</strong> Compare the native-ligand and mutant-ligand complexes. juProt highlights how the mutation alters the interaction network, which can help explain changes in binding affinity, drug efficacy, or resistance mechanisms. (e.g., comparing a wild-type kinase-inhibitor complex with a gatekeeper mutant-inhibitor complex).</p>

        <h3>2. Comparing Different Ligands to the Same Target</h3>
        <p><strong>Scenario:</strong> You have several drug candidates or chemical probes binding to the same protein target.</p>
        <p><strong>How juProt Helps:</strong> Compare Protein+LigandA with Protein+LigandB. juProt helps identify which ligand forms more/different interactions and which residues are key common or unique interaction partners, aiding in SAR studies and lead optimization.</p>

        <h3>3. Analyzing Ligand Binding to Different Protein Conformations or Isoforms</h3>
        <p><strong>Scenario:</strong> A protein exists in different states (e.g., active/inactive) or as different isoforms, and you have structures of a ligand bound to these variants.</p>
        <p><strong>How juProt Helps:</strong> Compare ProteinState1+Ligand with ProteinState2+Ligand. juProt can reveal how protein structural changes influence the interaction profile with a common ligand.</p>

        <h3>4. Validating Docking Poses</h3>
        <p><strong>Scenario:</strong> You have multiple potential binding poses for a ligand from molecular docking.</p>
        <p><strong>How juProt Helps:</strong> Compare the interaction profile of different docked poses or a docked pose against an experimental structure (if available) to assess consistency.</p>

        <h3>5. Educational Purposes</h3>
        <p><strong>Scenario:</strong> Teaching students about protein-ligand interactions.</p>
        <p><strong>How juProt Helps:</strong> Provides an easy-to-use tool for students to explore differential interactions without complex software or scripting.</p>

        <h3>Future Enhancements</h3>
        <p>juProt is an actively developed open-source project. Future versions aim to include:</p>
        <ul>
            <li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li>
            <li>Allowing analysis of more than two complexes.</li>
            <li>User-configurable parameters for interaction detection.</li>
        </ul>

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
            h1, h2 { color: #0056b3; }
            h1 { text-align: center; margin-bottom: 30px; }
            h2 { margin-top: 25px; border-bottom: 1px solid #eee; padding-bottom: 5px;}
            p, li { line-height: 1.6; }
            a { color: #007bff; text-decoration: none; }
            a:hover { text-decoration: underline; }
            .footer-link { text-align: center; margin-top: 30px; }
        </style>
    </head>
    <body>
        <div class="container">
        <div style="text-align: center; margin-bottom: 20px;"> 
            <img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height: 70px; margin-top: 10px;"> 
            </div>
            <h1>About juProt</h1>

            <p>juProt (<code>https://juprot.info/</code>) is an open-source web application designed to facilitate the rapid and user-friendly comparative analysis of protein-ligand interaction networks, covering the complete spectrum of non-covalent interactions including Hydrogen Bonds, Hydrophobic Contacts, Salt Bridges, Water Bridges, Pi-Stacking, and Halogen Bonds. Understanding how ligands interact with their protein targets, and how these interactions change due to mutations or when comparing different ligands, is fundamental in structural biology, bioinformatics, and drug discovery.</p>

            <h2>Motivation</h2>
            <p>While powerful tools exist for analyzing interactions in a single protein-ligand complex, comparing these interactions across two different complexes often requires manual data extraction, scripting, and collation of results. juProt aims to simplify this process, making comparative interaction analysis accessible to a broader range of researchers and students without requiring extensive computational expertise.</p>

            <h2>Core Technology</h2>
            <p>juProt is developed using the Julia programming language, leveraging the high-performance Genie.jl web framework for its backend and user interface. The core interaction detection is powered by the well-established Protein-Ligand Interaction Profiler (PLIP), a Python-based tool <a href="https://academic.oup.com/nar/article/43/W1/W443/2467865">(Sebastian et al., 2015)</a>. juProt interfaces with PLIP using the PythonCall.jl package, allowing seamless integration of PLIP's robust algorithms.</p>

            <h2>Current Features</h2>
            <ul>
                <li>Upload of two PDB files for comparison.</li>
                <li>Automated identification of potential ligands with user selection.</li>
                <li>Detection and quantification of full interactomes (H-bonds, Hydrophobic, Salt Bridges, etc.) for each complex.</li>
                <li>Generation of:
                    <ul>
                        <li>A comparative summary table (CSV) highlighting common and differential interactions and residues.</li>
                        <li>A detailed list of all interactions for both complexes (CSV).</li>
                        <li>A bar chart visually comparing interaction profiles per residue (PNG).</li>
                        <li>An on-page analytical summary of key findings.</li>
                    </ul>
                </li>
            </ul>

            <h2>Future Development</h2>
            <p>juProt is an ongoing project. Future development plans include:</p>
            <ul>
                <li>Enhanced visualization options, potentially including 2D comparative interaction diagrams.</li>
                <li>Allowing analysis of more than two complexes.</li>
                <li>User-configurable parameters for interaction detection.</li>
            </ul>

            <h2>Development Team</h2>
            <p>juProt was conceived and developed by:<br>
               <a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br>
               Deepak S P, MSc Biotechnology<br> 
               Siva V, MSc Biotechnology<br>
               Surya Sekaran, [PhD]
            </p>
            
            <p>We also acknowledge the developers of the core libraries used in juProt, including Julia, Genie.jl, PythonCall.jl, PLIP, OpenBabel, and CairoMakie.</p>

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