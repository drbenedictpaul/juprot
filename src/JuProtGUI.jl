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
     
        footer {
    text-align: center;
    margin-top: 40px;
    padding-top: 20px;
    border-top: 1px solid #eee;
    font-size: 0.9em;
    color: #6c757d;
}
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

        <footer>
    <p>Developed by<br>
    <a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br>
    <a href="https://deepakshankar810.github.io/portfolio/" target="_blank">Deepak S P</a><br>
    <a href="https://siva1106.github.io/website/" target="_blank">Siva V</a><br>
    <a href="https://www.linkedin.com/in/surya-s-09b655166/" target="_blank">Surya Sekaran, [PhD]
    </p>
</footer>

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
           
                footer {
    text-align: center;
    margin-top: 40px;
    padding-top: 20px;
    border-top: 1px solid #eee;
    font-size: 0.9em;
    color: #6c757d;
}
           
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

                <footer>
    <p>Developed by<br>
    <a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br>
    <a href="https://deepakshankar810.github.io/portfolio/" target="_blank">Deepak S P</a><br>
    <a href="https://siva1106.github.io/website/" target="_blank">Siva V</a><br>
    <a href="https://www.linkedin.com/in/surya-s-09b655166/" target="_blank">Surya Sekaran, [PhD]
    </p>
</footer>
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
        <!DOCTYPE html><html><head><title>Results</title><link rel="icon" href="/img/favicon.ico"><style>body{font-family:Arial,sans-serif;margin:40px;background:#f4f7f6;color:#333}.container{max-width:900px;margin:auto;background:#fff;padding:20px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,.1)}.btn{display:inline-block;padding:10px 15px;background:#007bff;color:#fff;border-radius:4px;text-decoration:none;margin-right:10px}.btn:hover{background:#0056b3}h1{text-align:center;color:#0056b3}h2{color:#007bff;border-bottom:1px solid #eee;padding-bottom:5px;margin-top:30px}img{max-width:100%;height:auto;display:block;margin:0 auto}pre{background:#e9ecef;padding:15px;border-radius:5px;white-space:pre-wrap}
        
        footer {
    text-align: center;
    margin-top: 40px;
    padding-top: 20px;
    border-top: 1px solid #eee;
    font-size: 0.9em;
    color: #6c757d;
}

        </style></head><body><div class="container">
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
        </div>
        
        <footer>
    <p>Developed by<br>
    <a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br>
    <a href="https://deepakshankar810.github.io/portfolio/" target="_blank">Deepak S P</a><br>
    <a href="https://siva1106.github.io/website/" target="_blank">Siva V</a><br>
    <a href="https://www.linkedin.com/in/surya-s-09b655166/" target="_blank">Surya Sekaran, [PhD]
    </p>
</footer>

        </body></html>
        """
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
        footer { text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #eee; font-size: 0.9em; color: #6c757d; }
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
            <li><b>Provide Input Complexes:</b> You have two options for each complex. Either upload PDB files directly from your computer, or enter the 4-character RCSB PDB IDs (e.g., <code>3EQM</code>) to fetch structures from the database.</li>
            <li><b>Load & Select Ligands:</b> Click "Load Ligands." On the next page, juProt will auto-detect potential ligands. Select the specific ligand of interest for each complex from the dropdown menu.</li>
            <li><b>Run Analysis & View Results:</b> Click "Run Analysis" to view the results page, which provides:
                <ul>
                    <li>An <b>Analytical Summary</b> for a quick overview of key statistics.</li>
                    <li>An <b>Interaction Count Profile</b> (bar chart) for a quantitative comparison of interactions per residue.</li>
                    <li>A <b>2D Spatial Pocket Projection</b>, which is a geometric map revealing the physical shape of the binding site.</li>
                    <li>Downloadable <b>CSV files</b> (Comparison Table and Detailed Interactions) for offline analysis.</li>
                </ul>
            </li>
        </ol>

        <h2>Applications & Use Cases for juProt</h2>
        <p>juProt is designed to provide rapid comparative insights into the <b>full non-covalent interactome</b> (Hydrogen bonds, Hydrophobic contacts, Salt Bridges, Pi-Stacking, etc.). Here are some scenarios where juProt can be particularly useful:</p>

        <h3>1. Understanding the Impact of Mutations</h3>
        <p><b>Scenario:</b> You have a wild-type protein-ligand structure and a mutant form (e.g., from a SNP or site-directed mutagenesis) bound to the same ligand.</p>
        <p><b>How juProt Helps:</b> Compare the native-ligand and mutant-ligand complexes. juProt highlights how the mutation alters the interaction network, which can help explain changes in binding affinity, drug efficacy, or resistance mechanisms.</p>

        <h3>2. Comparing Different Ligands to the Same Target</h3>
        <p><b>Scenario:</b> You have several drug candidates or chemical probes binding to the same protein target.</p>
        <p><b>How juProt Helps:</b> Compare Protein+LigandA with Protein+LigandB. juProt helps identify which ligand forms more or different interactions and which residues are key common or unique interaction partners, aiding in SAR studies and lead optimization.</p>

        <h3>3. Analyzing Ligand Binding to Different Protein Conformations or Isoforms</h3>
        <p><b>Scenario:</b> A protein exists in different states (e.g., active/inactive) or as different isoforms, and you have structures of a ligand bound to these variants.</p>
        <p><b>How juProt Helps:</b> Compare ProteinState1+Ligand with ProteinState2+Ligand. juProt can reveal how protein structural changes influence the interaction profile with a common ligand.</p>

        <h3>4. Validating Docking Poses</h3>
        <p><b>Scenario:</b> You have multiple potential binding poses for a ligand from molecular docking.</p>
        <p><b>How juProt Helps:</b> Compare the interaction profile of a docked pose against an experimental structure (if available) to assess interaction consistency.</p>

        <h3>5. Educational Purposes</h3>
        <p><b>Scenario:</b> Teaching students about protein-ligand interactions.</p>
        <p><b>How juProt Helps:</b> Provides an easy-to-use tool for students to explore differential interactions without complex software or scripting.</p>

        <div class="footer-link"><p><a href="/">Back to Home</a></p></div>
    </div>
    <footer>
        <p>Developed by<br>
        <a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br>
        <a href="https://deepakshankar810.github.io/portfolio/" target="_blank">Deepak S P</a>, MSc Biotechnology<br>
        <a href="https://siva1106.github.io/website/" target="_blank">Siva V</a>, MSc Biotechnology<br>
        Surya Sekaran, [PhD]
        </p>
    </footer>
</body>
</html>
""", forceparse=false)
   end
 
    route("/about") do
     html("""
     <!DOCTYPE html><html><head><title>About juProt</title><link rel="icon" href="/img/favicon.ico" type="image/x-icon"><style>body{font-family:Arial,sans-serif;margin:40px;background-color:#f4f7f6;color:#333}.container{max-width:800px;margin:auto;background-color:#fff;padding:20px 40px;border-radius:8px;box-shadow:0 0 10px rgba(0,0,0,0.1)}h1,h2{color:#0056b3}h1{text-align:center;margin-bottom:30px}h2{margin-top:25px;border-bottom:1px solid #eee;padding-bottom:5px}p,li{line-height:1.6}a{color:#007bff;text-decoration:none}a:hover{text-decoration:underline}.footer-link{text-align:center;margin-top:30px}blockquote{background:#e9ecef;border-left:5px solid #007bff;padding:15px;margin:20px 0;font-style:italic}footer{text-align:center;margin-top:40px;padding-top:20px;border-top:1px solid #eee;font-size:0.9em;color:#6c757d}</style></head><body><div class="container">
     <div style="text-align:center;margin-bottom:20px"><img src="/img/juProt_logo.png" alt="JuProt Logo" style="max-height:70px;margin-top:10px"></div>
     <h1>About juProt</h1>
     <p><b>juProt</b> is an open-source standalone application designed to facilitate the rapid and user-friendly comparative analysis of protein-ligand interaction networks, covering the complete spectrum of non-covalent interactions.</p>

     <h2>Motivation</h2>
     <p>While powerful tools exist for analyzing interactions in a single complex, comparing these interactions across two different complexes often requires manual scripting. juProt automates this process, making comparative interaction analysis accessible to a broader range of researchers.</p>
     
     <h2>Core Technology</h2>
     <p>juProt is developed using the <a href="https://julialang.org/" target="_blank">Julia programming language</a> and the <a href="https://genieframework.com/" target="_blank">Genie.jl</a> framework. The core interaction detection is powered by the well-established Protein-Ligand Interaction Profiler (<a href="https://doi.org/10.1093/nar/gkv315" target="_blank">PLIP</a>). The application is distributed as a self-contained Docker image.</p>

     <h2>Current Features</h2>
     <ul>
         <li>Upload of PDB files or direct fetch via PDB ID.</li>
         <li>Automated identification of ligands.</li>
         <li>Detection of the full non-covalent interactome (H-bonds, Hydrophobic, Salt Bridges, etc.).</li>
         <li>Generation of:
             <ul>
                 <li>A comparative summary table (CSV).</li>
                 <li>A detailed list of all interactions (CSV).</li>
                 <li>A <b>2D Spatial Pocket Projection</b> for geometric comparison.</li>
                 <li>A bar chart for quantitative comparison.</li>
                 <li>An on-page analytical summary.</li>
             </ul>
         </li>
     </ul>

     <h2>Development Team</h2>
     <p>juProt was conceived and developed by:<br>
     <a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br>
     <a href="https://deepakshankar810.github.io/portfolio/" target="_blank">Deepak S P</a>, MSc Biotechnology<br>
     <a href="https://siva1106.github.io/website/" target="_blank">Siva V</a>, MSc Biotechnology<br>
     Surya Sekaran, [PhD]</p>
     
     <p>We also acknowledge the developers of the core libraries used in juProt, including Julia, Genie.jl, PythonCall.jl, PLIP, OpenBabel, and Plots.jl.</p>

     <h2>Open Source & Citation</h2>
     <p>juProt is an open-source project. The source code is available on GitHub. If you use this tool in your research, please cite our publication:</p>
     <blockquote>
         <b>[Citation]</b><br>
         Benedict Christopher Paul, Deepak S P, Siva V, Surya Sekaran. juProt: A web application for comparative analysis of protein-ligand interactomes. In Silico Pharmacol. 2026 Feb 26;14(1):81. doi: 10.1007/s40203-026-00588-6. PMID: 41767850; PMCID: PMC12936254. <a href="https://doi.org/10.1007/s40203-026-00588-6" target="_blank">https://doi.org/10.1007/s40203-026-00588-6</a>
     </blockquote>

     <h2>Contact/Feedback</h2>
     <p>For questions, suggestions, or to report issues, please visit our GitHub Issues page: <a href="https://github.com/drbenedictpaul/juprot/issues" target="_blank">GitHub Issues</a>.</p>

     <div class="footer-link"><p><a href="/">Back to Home</a></p></div>
     </div>
     <footer><p>Developed by<br><a href="https://www.drpaul.cc/" target="_blank">Dr. Benedict Christopher Paul</a><br><a href="https://deepakshankar810.github.io/portfolio/" target="_blank">Deepak S P</a>, MSc Biotechnology<br><a href="https://siva1106.github.io/website/" target="_blank">Siva V</a>, MSc Biotechnology<br>Surya Sekaran, [PhD]</p></footer>
     </body></html>
     """, forceparse=false)
    end
 
 end # module JuProtGUI