# bootstrap.jl

cd(@__DIR__) # Ensure we are in the app's root directory
using Pkg
Pkg.activate(".") # Activate the project environment

# Load your main application module
# This assumes your main Genie application logic (routes, etc.) is in JuProtGUI.jl
# and JuProtGUI.jl is located inside the "src" folder.
include(joinpath("src", "JuProtGUI.jl")) 

using Genie # This makes Genie.config and Genie.up available

# --- Configure Genie for Cloud Run ---
# Directly modify the Genie.config object
Genie.config.server_host = "0.0.0.0"
Genie.config.server_port = parse(Int, get(ENV, "PORT", "8080")) # Default to 8080
Genie.config.run_as_server = true # Important for production
# Genie.config.websockets_port = parse(Int, get(ENV, "WS_PORT", "8081")) # If using websockets, uncomment and adjust

# --- Start Genie Server ---
println("Attempting to start Genie server with Genie.up on host: $(Genie.config.server_host), port: $(Genie.config.server_port)")

# Use Genie.up() to start the server.
# async=false is default when run_as_server=true, but explicitly setting it is fine.
Genie.up(Genie.config.server_port, Genie.config.server_host; async=false)

println("Genie server should be running if no errors occurred before this line.")