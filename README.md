# juProt: A Standalone Application for Comparative Protein-Ligand Analysis

<p align="center">
  <img src="/public/img/juProt_logo.png" alt="juProt Logo" width="200"/>
</p>

[![Version](https://img.shields.io/badge/version-v1.0-blue)](https://github.com/drbenedictpaul/juprot/releases/tag/v1.0)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**juProt** is an open-source, standalone application for the comparative analysis of protein-ligand interaction networks. It provides a user-friendly interface to compare two complexes, calculating and visualizing differences in **Hydrogen Bonds, Hydrophobic Contacts, Salt Bridges, Pi-Stacking, and more** using the PLIP engine.

**Status:** This project has been transitioned from a web service to a standalone Dockerized application to ensure its long-term availability. The core functionality is the same as described in our publication.

---

## Quick Start for Users

The primary way to use juProt is via the pre-packaged Docker application.

**➡️ Step 1: Download the application from the [Official Releases Page](https://github.com/drbenedictpaul/juprot/releases).**

**➡️ Step 2: Follow the installation guide below.**

---

### Installation and User Guide (for Linux)

This guide will walk you through setting up and running the juProt application.

#### Part 1: One-Time System Setup (Docker)

Before using juProt for the first time, your system needs one required piece of software: **Docker**. If you already have Docker installed, you can skip to Part 2.

1.  **Install Docker**
    *   **For Ubuntu/Debian:**
        ```bash
        sudo apt-get update && sudo apt-get install docker.io
        ```
    *   **For Fedora:**
        ```bash
        sudo dnf install docker && sudo systemctl start docker
        ```

2.  **Add Your User to the Docker Group** (Important!)
    This step allows you to run Docker commands without `sudo`.
    ```bash
    sudo usermod -aG docker ${USER}
    ```
    **CRITICAL:** You must **log out and log back in** for this change to take effect.

3.  **Verify Docker**
    After logging back in, run `docker run hello-world`. If you see a "Hello from Docker!" message, you are ready.

#### Part 2: Running juProt

1.  **Unpack the Application Package**
    In the folder where you downloaded the file, run:
    ```bash
    tar -xzvf juprot-linux-v1.0.tar.gz
    ```

2.  **Run juProt**
    To start the application, run the launcher script:
    ```bash
    ./run_juprot.sh
    ```
    The first time you run this, it will load the application into Docker. The script will then guide you.

3.  **Access the Application**
    Open your web browser and navigate to **[http://localhost:8888](http://localhost:8888)**.

4.  **Stopping the Application**
    To shut down the server, return to the terminal where it is running and press **Ctrl+C**. Use `./stop_juprot.sh` to clean up any background processes if needed.

---

### Key Features

*   **Comprehensive Interaction Analysis:** Detects a wide range of non-covalent interactions via the PLIP engine.
*   **Comparative Visualization:**
    *   **2D Spatial Pocket Map:** A geometric projection of the binding pocket.
    *   **Bar Charts:** Visual comparison of interaction counts per residue.
*   **Data Export:** Detailed CSV files for summary and in-depth analysis.
*   **Analytical Summary:** A text-based summary of key differences.

### Technology Stack

*   **Backend & Web Framework:** [Julia](https://julialang.org/) with [Genie.jl](https://genieframework.com/)
*   **Interaction Engine:** [PLIP](https://github.com/pharmai/plip) (Python-based)
*   **Distribution:** [Docker](https://www.docker.com/)

### Screenshots

![juProt Input Page](./public/img/juProt_input_page.png)
![juProt Select Ligands Page](./public/img/juProt_select_ligands_page.png)
![juProt Results Page](./public/img/aromatase_native_mutant.png)

---

### For Developers: Running from Source

If you wish to run the application directly from the source code for development:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/drbenedictpaul/juprot.git
    cd juprot
    ```

2.  **Run the application:**
    ```bash
    bash bin/server
    ```
    The first time you run this, Julia's package manager (`Pkg.jl`) and `CondaPkg` will automatically download and install all required Julia and Python dependencies. This may take several minutes. Subsequent startups will be fast.

### Citation

If you use juProt in your research, please cite our publication:

> *[Paul BC, S P D, V S, Sekaran S. juProt: A web application for comparative analysis of protein-ligand interactomes. In Silico Pharmacol. 2026 Feb 26;14(1):81. doi: 10.1007/s40203-026-00588-6. PMID: 41767850; PMCID: PMC12936254.]*

The source code for the version used in the publication is available under the `v1.0` tag.

### License

This project is licensed under the MIT License. See the `LICENSE` file for details.