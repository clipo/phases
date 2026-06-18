# Container for reproducing "Are the phases real? Distinguishing bounded
# interaction groups from spatially structured drift in central Mississippi
# Valley decorated ceramics."
#
# Build:  docker build -t phases .
# Run the full pipeline (writes into the container; mount to keep outputs):
#   docker run --rm -it \
#     -v "$PWD/output:/repo/output" -v "$PWD/figures:/repo/figures" \
#     phases ./run_all.sh
# Interactive shell:
#   docker run --rm -it phases bash
#
# The conda environment is created from environment.yml in its own layer, so it
# is cached and only rebuilt when environment.yml changes. For a bit-for-bit
# environment, pin the base image to a digest (FROM condaforge/miniforge3@sha256:...)
# and create the env from conda-lock.yml instead of environment.yml.
FROM condaforge/miniforge3:latest

WORKDIR /repo

# 1) Build the conda environment (cached unless environment.yml changes).
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# 2) Copy the project and install the package into the environment.
COPY . .
RUN conda run -n phases pip install -e .

# Put the environment on PATH so `python`, `pandoc`, and `./run_all.sh` use it.
ENV PATH=/opt/conda/envs/phases/bin:$PATH
ENV PYTHONUNBUFFERED=1

CMD ["bash"]
