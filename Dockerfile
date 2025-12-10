# Stage 1: The Builder - Compiles the application in a clean environment
FROM julia:1.10 AS builder

ENV JULIA_PROJECT=/app
ENV GENIE_ENV=prod
WORKDIR /app

# Copy ONLY Project.toml and CondaPkg.toml. DO NOT copy Manifest.toml.
COPY Project.toml CondaPkg.toml ./

# This will create a fresh, Debian-compatible Manifest.toml inside the container
RUN julia -e 'using Pkg; Pkg.instantiate()'

# Now copy the rest of the source code
COPY . .

# Pre-compile everything
RUN julia --project -e 'using Pkg; Pkg.precompile()'


# Stage 2: The Runner - A minimal image to run the app
FROM julia:1.10

ENV JULIA_PROJECT=/app
ENV GENIE_ENV=prod
ENV GKSwstype=100
WORKDIR /app

# Copy the entire pre-built application from the builder stage
COPY --from=builder /app .
# Copy the pre-compiled packages and the Conda environment
COPY --from=builder /usr/local/julia/packages /usr/local/julia/packages
COPY --from=builder /app/.CondaPkg /app/.CondaPkg

# The command to run the server
CMD ["julia", "--project=.", "bootstrap.jl"]
