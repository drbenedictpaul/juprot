# --- STAGE 1: Python Environment Builder ---
# We use an official micromamba image which is guaranteed to work.
# We name this stage "conda_builder" so we can refer to it later.
FROM mambaorg/micromamba:latest as conda_builder

# Define the installation path for our environment
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Create the 'bio' environment with all Python dependencies.
# We use --prefix to install it to a specific, predictable location.
RUN micromamba create -y --prefix /opt/conda/envs/bio -c conda-forge \
    python=3.11 \
    openbabel \
    plip \
    numpy


# --- STAGE 2: Julia Application Builder ---
# This is the final image that will be distributed.
FROM julia:1.12

WORKDIR /app

# Install only the necessary system libraries for Julia and Plots.jl
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libgl1 \
    libxrender1 \
    libxext6 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Copy the entire pre-built Python environment from the first stage.
COPY --from=conda_builder /opt/conda/envs/bio /opt/conda/envs/bio

# --- Set all necessary environment variables ---
# Tell Julia to use the Python we just copied.
ENV JULIA_CONDAPKG_BACKEND=Null
ENV JULIA_PYTHONCALL_EXE=/opt/conda/envs/bio/bin/python

# Headless Plotting Mode
ENV GKSwstype=100

# Set the Genie environment to production
ENV GENIE_ENV=prod

# --- THE KEY FIX ---
# Set the PORT environment variable that your bootstrap.jl is looking for.
ENV PORT=8888

# --- JULIA SETUP ---
# Copy project files first to leverage Docker's layer caching
COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.instantiate()'

# Copy the rest of the application code
COPY . .

# Precompile the full application. This will now find the Python packages.
RUN julia --project=. -e 'using Pkg; Pkg.precompile()'

# --- RUNTIME CONFIGURATION ---
# Expose the port that the application will run on.
EXPOSE 8888

# The bootstrap.jl file will now automatically pick up the PORT=8888 variable.
CMD ["julia", "--project=.", "bootstrap.jl"]