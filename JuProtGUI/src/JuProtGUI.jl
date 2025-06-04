# JuProtGUI/src/JuProtGUI.jl
module JuProtGUI
using Genie, Genie.Router, Genie.Renderer, Genie.Renderer.Html, Genie.Requests
using BioStructures, CSV, DataFrames, Plots, LinearAlgebra
using Genie.Assets

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
    native_path = joinpath(output_dir, "native.cif")
    mutated_path = joinpath(output_dir, "mutated.cif")
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
            for p_atom in protein_atoms
                for l_atom in ligand_atoms
                    dist = hbond_distance(p_atom, l_atom)
                    if dist < threshold
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
            nonbonded = detect_nonbonded(protein_atoms, ligand_atoms)
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
        <title>juProt GUI</title>
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
                    <label for="native_cif">Native CIF File:</label>
                    <input type="file" id="native_cif" name="native_cif" accept=".cif" required>
                </div>
                <div class="form-group">
                    <label for="native_ligand">Native Ligand Residue Name (e.g., ASD):</label>
                    <input type="text" id="native_ligand" name="native_ligand" required>
                </div>
                <div class="form-group">
                    <label for="mutated_cif">Mutated CIF File:</label>
                    <input type="file" id="mutated_cif" name="mutated_cif" accept=".cif" required>
                </div>
                <div class="form-group">
                    <label for="mutated_ligand">Mutated Ligand Residue Name (e.g., TES):</label>
                    <input type="text" id="mutated_ligand" name="mutated_ligand" required>
                </div>
                <button type="submit">Run Analysis</button>
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

end