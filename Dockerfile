# --- STAGE 1: Python Environment Builder ---
FROM mambaorg/micromamba:latest as conda_builder
ENV MAMBA_ROOT_PREFIX=/opt/conda
RUN micromamba create -y --prefix /opt/conda/envs/bio -c conda-forge \
    python=3.11 \
    openbabel \
    plip \
    numpy

# --- STAGE 2: Julia Application Builder ---
FROM julia:1.12
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libgl1 \
    libxrender1 \
    libxext6 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Copy the pre-built Python environment
COPY --from=conda_builder /opt/conda/envs/bio /opt/conda/envs/bio

# Set environment variables
ENV JULIA_CONDAPKG_BACKEND=Null
ENV JULIA_PYTHONCALL_EXE=/opt/conda/envs/bio/bin/python
ENV GKSwstype=100
ENV GENIE_ENV=prod
ENV PORT=8888

# --- CORRECTED AND ROBUST JULIA SETUP ---
# 1. Copy project files for dependency installation
COPY Project.toml Manifest.toml ./

# 2. Install all Julia dependencies. This step will be cached.
RUN julia -e 'using Pkg; Pkg.instantiate()'

# 3. Copy THE ENTIRE application context (all source code, public files, etc.)
COPY . .

# 4. Precompile the FULL application now that all files are present.
# This is the slow step, but its cache will only be broken if source code changes.
RUN julia --project=. -e 'using Pkg; Pkg.precompile()'

# --- RUNTIME CONFIGURATION ---
EXPOSE 8888
CMD ["julia", "--project=.", "bootstrap.jl"]