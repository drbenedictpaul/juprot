FROM julia:1.11.0

WORKDIR /app
COPY . /app

RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

EXPOSE 8000
CMD ["julia", "--project=.", "-e", "using Genie; Genie.loadapp(); Genie.up(8000, \"0.0.0.0\", async=false)"]
