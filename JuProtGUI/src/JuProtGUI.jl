# JuProtGUI/src/JuProtGUI.jl
module JuProtGUI
using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
using BioStructures, CSV, DataFrames, Plots, LinearAlgebra
using Genie.Assets

const SETTINGS = Dict{Symbol, Any}()

include("../lib/utils.jl")
include("../lib/output_results.jl")
include("../lib/detect_hbonds.jl")
include("../lib/detect_nonbonded.jl")
include("../lib/detect_pi_pi.jl")

function process_cif_files(native_file, mutated_file, native_ligand, mutated_ligand)
    complex_results = Dict{String, Dict}()
    structures = Dict{String, Any}()
    output_dir = "public/outputs"
    mkpath(output_dir)
    native_path = joinpath(output_dir, "first_complex.cif")
    mutated_path = joinpath(output_dir, "second_complex.cif")
    write(native_path, native_file.data)
    write(mutated_path, mutated_file.data)
    cif_files = [native_path, mutated_path]
    ligand_resnames = [native_ligand, mutated_ligand]

    for (i, cif_file) in enumerate(cif_files)
        try
            structure = read(cif_file, MMCIFFormat)
            models = collectmodels(structure)
            model = structure[1]
            protein_atoms = collectatoms(model, standardselector)
            non_protein_atoms = collectatoms(model, atom -> !standardselector(atom) && resname(atom) != "HOH" && resname(atom) != "WAT")
            unique_resnames = unique(resname(atom) for atom in non_protein_atoms)
            ligand_resname = ligand_resnames[i]
            if !(ligand_resname in unique_resnames)
                return Dict("error" => "Ligand '$ligand_resname' not found in $cif_file")
            end
            ligand_selector(atom) = resname(atom) == ligand_resname
            ligand_atoms = collectatoms(model, ligand_selector)
            close_contacts = []
            threshold = 4.0
            hbond_range = get(SETTINGS, :hbond_distance_range, (2.0, 3.5))
            nonbonded_range = get(SETTINGS, :nonbonded_distance_range, (3.0, 4.0))
            for p_atom in protein_atoms
            for l_atom in ligand_atoms
                dist = hbond_distance(p_atom, l_atom)
                    if hbond_range[1] <= dist <= hbond_range[2]
                        push!(close_contacts, (p_atom, l_atom, dist, "Close Contact"))
                    end
                end
            end
            interacting_residues = Dict{String, Int}()
            for (p_atom, _, _, _) in close_contacts
                res_key = "$(resname(p_atom)) $(resnumber(p_atom))"
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
            hbonds = detect_hbonds(protein_atoms, ligand_atoms)
            nonbonded = detect_nonbonded(protein_atoms, ligand_atoms, nonbonded_range[1], nonbonded_range[2])
            pi_pi = detect_pi_pi(protein_atoms, ligand_atoms)
            complex_name = basename(cif_file)
            complex_results[complex_name] = Dict(
                :close_contacts => close_contacts,
                :hbonds => hbonds,
                :nonbonded => nonbonded,
                :pi_pi => pi_pi,
                :interacting_residues => interacting_residues,
                :max_residue => max_residue,
                :max_count => max_count,
                :ligand_resname => ligand_resname
            )
            structures[complex_name] = structure
        catch e
            return Dict("error" => "Failed to process $cif_file: $(sprint(showerror, e))")
        end
    end

    output_file = joinpath(output_dir, "comparison_table.csv")
    plot_file = joinpath(output_dir, "residue_interactions.png")
    try
        compare_and_save_results(complex_results, output_file, plot_file, 
                                structures[basename(cif_files[1])], structures[basename(cif_files[2])])
        summary = capture_analytical_summary(complex_results)
        return Dict(
            "comparison_table" => output_file,
            "detailed_interactions" => joinpath(output_dir, "detailed_interactions.csv"),
            "plot_file" => plot_file,
            "summary" => summary
        )
    catch e
        return Dict("error" => "Analysis failed: $(sprint(showerror, e))")
    end
end

function capture_analytical_summary(complex_results)
    output = Pipe()
    redirect_stdout(output) do
        print_analytical_summary(complex_results)
    end
    close(output.in)
    return String(read(output))
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
            input[type=file], input[type=text] { width: 100%; padding: 8px; }
            button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; }
            button:hover { background: #0056b3; }
            pre { background: #f8f9fa; padding: 10px; border-radius: 5px; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>juProt: Protein-Ligand Interaction Analysis</h1>
            <form action="/analyze" method="post" enctype="multipart/form-data">
                <div class="form-group">
                    <label for="native_cif">Native or First Complex CIF File:</label>
                    <input type="file" id="native_cif" name="native_cif" accept=".cif" required>
                </div>
                <div class="form-group">
                    <label for="native_ligand">Native or First Complex Ligand Name (e.g., ASD):</label>
                    <input type="text" id="native_ligand" name="native_ligand" required>
                </div>
                <div class="form-group">
                    <label for="mutated_cif">Mutated or Second Complex CIF File:</label>
                    <input type="file" id="mutated_cif" name="mutated_cif" accept=".cif" required>
                </div>
                <div class="form-group">
                    <label for="mutated_ligand">Mutated or Second Complex Ligand Name (e.g., TES):</label>
                    <input type="text" id="mutated_ligand" name="mutated_ligand" required>
                </div>
                <button type="submit">Run Analysis</button>
                <p><a href="/how-to-use">How to Use</a></p>
                <p><a href="/settings">Settings</a></p>
            </form>
        </div>
    </body>
    </html>
    """)
end

route("/analyze", method=POST) do
    native_file = filespayload("native_cif")
    mutated_file = filespayload("mutated_cif")
    native_ligand = postpayload(:native_ligand)
    mutated_ligand = postpayload(:mutated_ligand)
    result = process_cif_files(native_file, mutated_file, native_ligand, mutated_ligand)
    if haskey(result, "error")
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
            <p>Comparison of interactions for the First Complex and Second Complex:</p>
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
            <p>juProt is a web-based tool for comparing protein-ligand interactions between two complexes, including native (e.g., from X-ray crystallography or NMR) or docked complexes (e.g., from molecular docking of native or mutated proteins).</p>
            <h2>Protocol</h2>
            <ol>
                <li><strong>Prepare CIF Files</strong>: Obtain two protein-ligand complex CIF files (e.g., from PDB or docking software).</li>
                <li><strong>Access the App</strong>: Visit <a href="/">juProt</a>.</li>
                <li><strong>Upload Files</strong>:
                    <ul>
                        <li>Upload the first CIF file under "Native or First Complex CIF File".</li>
                        <li>Enter the ligand residue name (e.g., ASD) for the first complex.</li>
                        <li>Upload the second CIF file under "Mutated or Second Complex CIF File".</li>
                        <li>Enter the ligand residue name (e.g., TES) for the second complex.</li>
                    </ul>
                </li>
                <li><strong>Run Analysis</strong>: Click "Run Analysis" to compare interactions (hydrogen bonds, non-bonded contacts, pi-pi interactions).</li>
                <li><strong>View Results</strong>:
                    <ul>
                        <li>Download the comparison table (CSV) for interaction counts.</li>
                        <li>Download detailed interactions (CSV) for specific bonds.</li>
                        <li>View or download the bar chart (PNG) comparing residue interactions.</li>
                        <li>Read the analytical summary for insights.</li>
                    </ul>
                </li>
            </ol>
            <h2>Applications</h2>
            <ul>
                <li><strong>Structural Biology</strong>: Compare native and mutated protein-ligand complexes to study mutation effects.</li>
                <li><strong>Drug Design</strong>: Analyze docked complexes to evaluate ligand binding differences.</li>
                <li><strong>Protein Engineering</strong>: Assess interaction changes in engineered proteins.</li>
            </ul>
            <p><a href="/">Back to Home</a></p>
        </div>
    </body>
    </html>
    """)
end

route("/settings") do
    html("""
    <!DOCTYPE html>
    <html>
    <head>
        <title>juProt: Settings</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .container { max-width: 800px; margin: auto; }
            h1, h2 { color: #007bff; }
            .form-group { margin-bottom: 20px; }
            label { display: block; margin-bottom: 5px; }
            input[type=number] { width: 100%; padding: 8px; }
            button { padding: 10px 20px; background: #007bff; color: white; border: none; cursor: pointer; }
            button:hover { background: #0056b3; }
            a { color: #007bff; text-decoration: none; }
            a:hover { text-decoration: underline; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>juProt: Settings</h1>
            <p>Customize the bond distance ranges for interaction analysis.</p>
            <form action="/settings" method="post">
                <div class="form-group">
                    <label for="hbond_distance">H-bond Distance Range (Å):</label>
                    <input type="number" id="hbond_distance_min" name="hbond_distance_min" step="0.1" min="0" max="5" value="2.0" required> to 
                    <input type="number" id="hbond_distance_max" name="hbond_distance_max" step="0.1" min="0" max="5" value="3.5" required>
                </div>
                <div class="form-group">
                    <label for="nonbonded_distance">Non-bonded Distance Range (Å):</label>
                    <input type="number" id="nonbonded_distance_min" name="nonbonded_distance_min" step="0.1" min="0" max="5" value="3.0" required> to 
                    <input type="number" id="nonbonded_distance_max" name="nonbonded_distance_max" step="0.1" min="0" max="5" value="4.0" required>
                </div>
                <button type="submit">Save Settings</button>
            </form>
            <p><a href="/">Back to Home</a></p>
        </div>
    </body>
    </html>
    """)
end

route("/settings", method=POST) do
    hbond_min = parse(Float64, postpayload(:hbond_distance_min, "2.0"))
    hbond_max = parse(Float64, postpayload(:hbond_distance_max, "3.5"))
    nonbonded_min = parse(Float64, postpayload(:nonbonded_distance_min, "3.0"))
    nonbonded_max = parse(Float64, postpayload(:nonbonded_distance_max, "4.0"))
    SETTINGS[:hbond_distance_range] = (hbond_min, hbond_max)
    SETTINGS[:nonbonded_distance_range] = (nonbonded_min, nonbonded_max)
    html("<h1>Settings Saved</h1><p>Bond distance ranges updated successfully. <a href='/'>Return to Home</a></p>")
end

end
