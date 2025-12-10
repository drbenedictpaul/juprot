cd(@__DIR__)

using Pkg
Pkg.activate(".")

# Explicitly set the environment to Production
ENV["GENIE_ENV"] = "prod"

using Genie
# This command loads the app, including the correct prod.jl config
Genie.loadapp()

# This command starts the server using the loaded config
Genie.startup()