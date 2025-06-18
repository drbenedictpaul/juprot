# Use an official Julia image as a parent image
# Or your preferred stable Julia 1.x version
FROM julia:1.11.5

# Set the working directory in the container
WORKDIR /app

# Install system dependencies for Python, OpenBabel, and PLIP
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    cmake \         
    build-essential \
    swig \
    pkg-config \
    libopenbabel-dev \
    && rm -rf /var/lib/apt/lists/*
    

# Create a Python virtual environment
ENV VENV_PATH=/opt/venv
RUN python3 -m venv $VENV_PATH 

# Add venv to PATH for subsequent RUN commands
ENV PATH="$VENV_PATH/bin:$PATH"

# Install PLIP and OpenBabel Python bindings into the venv
RUN pip install --no-cache-dir plip openbabel-wheel

# --- Julia Application Setup ---
COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'
COPY . .

# --- Environment Variables for Runtime ---
ENV JULIA_PYTHONCALL_EXE=$VENV_PATH/bin/python
ENV JULIA_CONDAPKG_BACKEND=Null
ENV PYTHON=""
ENV JULIA_PYTHONCALL_LIB=/usr/lib/x86_64-linux-gnu/libpython3.11.so

EXPOSE 8080
CMD ["julia", "--project=.", "bootstrap.jl"]