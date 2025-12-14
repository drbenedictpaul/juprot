# 1. Base Image
FROM julia:1.12

WORKDIR /app

# 2. Install System Basics + GRAPHICS LIBRARIES (Required for Plots.jl)
# Added: libgl1, libxrender1, libxext6, libglib2.0-0
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    bzip2 \
    ca-certificates \
    git \
    libgl1 \
    libxrender1 \
    libxext6 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# 3. Install Micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# 4. Create the Bioinformatics Environment
ENV MAMBA_ROOT_PREFIX=/opt/conda
RUN ./bin/micromamba create -y -n bio -c conda-forge \
    python=3.11 \
    openbabel \
    plip \
    numpy

# 5. Tell Julia to use this specific Python
ENV JULIA_CONDAPKG_BACKEND=Null
ENV JULIA_PYTHONCALL_EXE=/opt/conda/envs/bio/bin/python

# Headless Plotting Mode
ENV GKSwstype=100

# 6. Julia Setup
COPY Project.toml ./
RUN julia -e 'using Pkg; Pkg.Registry.add("General"); Pkg.activate("."); Pkg.resolve(); Pkg.instantiate()'

# 7. Copy Application Code
COPY . .

# 8. Precompile
RUN julia --project=. -e 'using Pkg; Pkg.precompile()'

# 9. Runtime Configuration
ENV GENIE_ENV=prod
EXPOSE 8080
CMD ["julia", "--project=.", "bootstrap.jl"]