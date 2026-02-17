FROM dolfinx/dolfinx:v0.7.2

# Install additional dependencies
RUN pip install --no-cache-dir \
    gmsh \
    pyvista \
    meshio \
    h5py \
    tqdm \
    pydantic \
    pyyaml \
    black \
    isort \
    mypy \
    flake8 \
    pytest \
    pytest-xdist \
    jupyter \
    ipywidgets \
    matplotlib

# Set up working directory
WORKDIR /workspace

# Copy source code
COPY . /workspace/

# Install package in editable mode
RUN pip install -e .

# Expose Jupyter port
EXPOSE 8888

# Default command
CMD ["/bin/bash"]
