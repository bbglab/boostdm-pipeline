# BoostDM

BoostDM is a machine-learning framework designed for predicting driver mutations in cancer. It leverages features like mutational clustering, functional domain enrichment, and nucleotide conservation to identify mutations with oncogenic potential.

---

## Getting Started

This guide will help you set up your development environment using `uv` and install the necessary dependencies to start working on the project.

### Prerequisites

Before you begin, ensure the following software is installed on your system -read the [HOWTO](/containers_build/HOW_TO.txt) for more info:

- **Python 3.9.20 or later**
- **uv** (Universal Virtual environment package manager)  
  Follow the [official uv installation guide](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer) to set it up.
- **A C/C++ compiler** This will be needed to install boostDM dependancies

### Installation

1. Clone the Repository:
   ```bash
   git clone https://github.com/bbglab/boostdm-pipeline.git
   cd containers_build
   ```

2. Install Project Dependancies and package:
   ```bash
   uv sync
   ```
   This will set up a `.venv/` directory to manage your project-specific environment.

### Development Workflow

1. **Adding Dependencies**:
   To add a new package to the project, use:
   ```bash
   uv add package-name
   ```
   Example:
   ```bash
   uv add numpy
   ```
   This updates the `pyproject.toml` and locks the dependency versions in `uv.lock`.

2. **Running BoostDM scripts via cmd line**:
   To run scripts or commands within the `uv` environment, ensure it is activated. For example:
   ```bash
   uv run boostdm --help
   ```

3. **Deactivating the Environment**:
   When you're done working, deactivate the environment:
   ```bash
   uv deactivate
   ```

