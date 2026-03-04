# juProt: A Standalone Application for Comparative Protein-Ligand Analysis

<p align="center">
  <img src="/public/img/juProt_logo.png" alt="juProt Logo" width="200"/>
</p>

[![Version](https://img.shields.io/badge/version-v1.0-blue)](https://github.com/drbenedictpaul/juprot/releases/tag/v1.0)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**juProt** is an open-source, standalone application for the comparative analysis of protein-ligand interaction networks. It provides a user-friendly interface to compare two complexes, calculating and visualizing differences in **Hydrogen Bonds, Hydrophobic Contacts, Salt Bridges, Pi-Stacking, and more** using the PLIP engine.

---

### **Status Update: Transition to a Standalone Application**

This project was originally published as a web service. To ensure its **long-term availability** and provide a more robust, cost-effective distribution model, **juProt has been transitioned into a standalone, Dockerized application** that runs locally on your machine.

The core scientific functionality and user interface remain the same as described in our publication, with the significant upgrade of supporting the full non-covalent interactome.

---

## Quick Start for Users (Linux & Windows)

The primary way to use juProt is via the pre-packaged Docker application.

**➡️ Step 1: Download the application from the [Official Releases Page](https://github.com/drbenedictpaul/juprot/releases).**

**➡️ Step 2: Follow the detailed installation guide below.**

---

## Installation and User Guide

This guide will walk you through setting up and running the juProt application.

### **Instructions for Linux Users**

#### Part 1: One-Time System Setup (Docker)
If you already have Docker installed, you can skip to Part 2.

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
    ```bash
    sudo usermod -aG docker ${USER}
    ```
    **CRITICAL:** You must **log out and log back in** for this change to take effect.

3.  **Verify Docker**
    After logging back in, run `docker run hello-world`. If you see a "Hello from Docker!" message, you are ready.

#### Part 2: Running juProt
1.  **Unpack the Application Package**
    In the folder where you downloaded the file, run: `tar -xzvf juprot-linux-v1.0.tar.gz`

2.  **Run juProt**
    `./run_juprot.sh`

3.  **Access the Application**
    Open your web browser and navigate to **[http://localhost:8888](http://localhost:8888)**.

### **Instructions for Windows Users (via WSL2)**

1.  **Install WSL & Docker Desktop**
    *   First, install the Windows Subsystem for Linux (WSL) by running `wsl --install` in an Administrator PowerShell. Reboot when prompted.
    *   Next, download and install **[Docker Desktop for Windows](https://www.docker.com/products/docker-desktop)**, ensuring the option "Use WSL 2" is selected during setup.

2.  **Open Your Linux Terminal**
    *   From your Start Menu, open the **Ubuntu** application. You are now inside a Linux terminal.

3.  **Navigate to your Downloads**
    *   Your Windows files are accessible at `/mnt/c/`. To get to your Downloads folder, run:
        ```bash
        cd /mnt/c/Users/YourUsername/Downloads
        ```
        *(Replace `YourUsername` with your actual Windows username).*

4.  **Unpack and Run**
    *   The steps are now identical to Linux.
    *   `tar -xzvf juprot-linux-v1.0.tar.gz`
    *   `./run_juprot.sh`

5.  **Access the Application**
    *   Open your regular Windows browser (Chrome, Edge, etc.) and go to **[http://localhost:8888](http://localhost:8888)**.

---

### Key Features & Technology

*   **Comprehensive Interaction Analysis:** Powered by the [PLIP](https://github.com/pharmai/plip) engine.
*   **Comparative Visualization:**
    *   **2D Spatial Pocket Map:** A geometric projection of the binding pocket.
    *   **Bar Charts:** Visual comparison of interaction counts per residue.
*   **Backend:** [Julia](https://julialang.org/) with [Genie.jl](https://genieframework.com/).
*   **Distribution:** [Docker](https://www.docker.com/).

### Screenshots

![juProt Input Page](./public/img/juProt_input_page.png)
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
    The first time you run this, all Julia and Python dependencies will be automatically installed, which may take several minutes.

### Citation

If you use juProt in your research, please cite our publication:

> *[Paul BC, S P D, V S, Sekaran S. juProt: A web application for comparative analysis of protein-ligand interactomes. In Silico Pharmacol. 2026 Feb 26;14(1):81. doi: 10.1007/s40203-026-00588-6. PMID: 41767850; PMCID: PMC12936254.]*

The source code for the version used in the publication is available under the `v1.0` tag.

### License

This project is licensed under the MIT License. See the `LICENSE` file for details.