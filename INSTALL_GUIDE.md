# juProt v1.0 - Installation and User Guide

This guide will walk you through the setup and execution of the juProt application on either a native Linux computer or on Windows via the Windows Subsystem for Linux (WSL).

---

## Instructions for Linux Users

#### Part 1: One-Time System Setup (Docker)

If you already have Docker installed, you can skip to Part 2.

1.  **Install Docker**
    *   **For Ubuntu/Debian:** `sudo apt-get update && sudo apt-get install docker.io`
    *   **For Fedora:** `sudo dnf install docker && sudo systemctl start docker`

2.  **Add Your User to the Docker Group** (Important!)
    ```bash
    sudo usermod -aG docker ${USER}
    ```
    **CRITICAL:** You must **log out and log back in** for this change to take effect.

3.  **Verify Docker**
    After logging back in, run `docker run hello-world`. If you see a "Hello from Docker!" message, you are ready.

#### Part 2: Running juProt

1.  **Run juProt**
    To start the application, simply run the launcher script from this folder:
    ```bash
    ./run_juprot.sh
    ```
    The first time you run this, it will load the application into Docker. The script will then guide you.

2.  **Access the Application**
    Open your web browser and navigate to **[http://localhost:8888](http://localhost:8888)**.

---

## Instructions for Windows Users (via WSL2)

1.  **Install WSL & Docker Desktop**
    *   First, install the Windows Subsystem for Linux (WSL) by running `wsl --install` in an Administrator PowerShell. Reboot when prompted.
    *   Next, download and install **[Docker Desktop for Windows](https://www.docker.com/products/docker-desktop)**, ensuring the option "Use WSL 2" is selected during setup.

2.  **Open Your Linux Terminal**
    *   From your Start Menu, open the **Ubuntu** application. You are now inside a Linux terminal.

3.  **Navigate to this Folder**
    *   Your Windows files are accessible at `/mnt/c/`. To get to your Downloads folder, for example, run:
        ```bash
        cd /mnt/c/Users/YourUsername/Downloads/juprot-folder
        ```
        *(Replace with the actual path to this folder).*

4.  **Run juProt**
    *   The steps are now identical to Linux: `./run_juprot.sh`

5.  **Access the Application**
    *   Open your regular Windows browser and go to **[http://localhost:8888](http://localhost:8888)**.

---

### Stopping the Application

To shut down the server, return to the terminal where it is running and press **Ctrl+C**.

If the container is ever left running in the background, you can use the failsafe script to stop it: `./stop_juprot.sh`
