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

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libgl1 \
    libxrender1 \
    libxext6 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=conda_builder /opt/conda/envs/bio /opt/conda/envs/bio

ENV JULIA_CONDAPKG_BACKEND=Null
ENV JULIA_PYTHONCALL_EXE=/opt/conda/envs/bio/bin/python
ENV GKSwstype=100

COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.instantiate()'

COPY . .
RUN julia --project=. -e 'using Pkg; Pkg.precompile()'

ENV GENIE_ENV=prod
EXPOSE 8888

# --- CORRECTED COMMAND ---
# We use 'bash' to execute the server script, which then starts Julia.
CMD ["bash", "bin/server", "0.0.0.0", "8888"]
