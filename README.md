# juProt: A Standalone Application for Comparative Protein-Ligand Analysis

<p align="center">
  <img src="/public/img/juProt_logo.png" alt="juProt Logo" width="200"/>
</p>

[![Version](https://img.shields.io/badge/version-v1.0-blue)](https://github.com/drbenedictpaul/juprot/releases/tag/v1.0)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**juProt** is an open-source, standalone application designed to streamline the comparative analysis of protein-ligand interaction networks. It helps researchers quickly identify differences and commonalities in how ligands bind to proteins by analyzing **Hydrogen Bonds, Hydrophobic Contacts, Salt Bridges, Pi-Stacking, and more.**

---

### **Status Update: Transition to a Standalone Application**

This project was originally published as a web service. To ensure its **long-term availability** and provide a more robust, cost-effective distribution model, **juProt has been transitioned into a standalone, Dockerized application** that runs locally on your machine.

---

## Quick Start for Users (Linux & Windows)

**➡️ Step 1: Download the application from the [Official Releases Page](https://github.com/drbenedictpaul/juprot/releases).**

**➡️ Step 2: The downloaded package (`.tar.gz`) contains a file named `INSTALL_GUIDE.md`. Follow the simple instructions in that guide to run the application.**

---

### Key Features

*   **Comprehensive Interaction Analysis:** Powered by the [PLIP](https://github.com/pharmai/plip) engine.
*   **Comparative Visualization:**
    *   **2D Spatial Pocket Map:** A geometric projection of the binding pocket.
    *   **Bar Charts:** Visual comparison of interaction counts per residue.
*   **Data Export:** Detailed CSV files for summary and in-depth analysis.

### For Developers: Running from Source

If you wish to run the application directly from the source code for development:

1.  **Clone the repository:** `git clone https://github.com/drbenedictpaul/juprot.git`
2.  **Run the application:** `cd juprot && bash bin/server`

The first time you run this, all Julia and Python dependencies will be automatically installed.

### Citation

If you use juProt in your research, please cite our publication:

> *[Your Full Citation Here, e.g., Paul, B.C., et al. juProt: a web application... in silico Pharmacology (2024). DOI: XXXXX]*

### License

This project is licensed under the MIT License. See the `LICENSE` file for details.
