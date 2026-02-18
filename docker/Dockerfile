FROM dolfinx/dolfinx:v0.7.2

# Install h5py from source to match system HDF5
# This prevents version mismatch errors
RUN pip install --no-cache-dir --no-binary=h5py h5py

# Install additional dependencies
RUN pip install --no-cache-dir \
    gmsh \
    pyvista \
    meshio \
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
