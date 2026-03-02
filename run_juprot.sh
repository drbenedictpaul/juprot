#!/bin/bash

# --- PART 1: Banner & Info ---
echo "========================================================================"
echo "          juProt v1.0 - Comparative Protein-Ligand Analyzer"
echo
echo "     Published in: in silico Pharmacology (Springer Nature)"
echo "     Developed by: Dr. Benedict Christopher Paul & Team"
echo "========================================================================"
echo

# --- PART 2: System Checks ---
IMAGE_NAME="juprot:latest"
CONTAINER_NAME="juprot-container"

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
  echo "❌ Error: Docker does not seem to be running. Please start Docker and try again."
  exit 1
fi

# Check if the user has loaded the application image from the .tar file
if ! docker image inspect $IMAGE_NAME > /dev/null 2>&1; then
    echo "🔎 juProt image not found. Loading it from juprot.tar..."
    echo "   This is a one-time setup and may take a minute."
    
    if [ -f "juprot.tar" ]; then
        docker load < juprot.tar
        echo "✅ Image loaded successfully."
    else
        echo "❌ Error: juprot.tar file not found. Make sure it is in the same folder as this script."
        exit 1
    fi
    echo
fi

# Check if the container is already running from a previous session
if [ "$(docker ps -q -f name=$CONTAINER_NAME)" ]; then
    echo "🟢 juProt is already running in another terminal."
    echo "   Access it at http://localhost:8888"
    echo "   To stop it, run ./stop_juprot.sh"
    exit 0
fi

# --- PART 3: LAUNCH ---
read -p "Press [Enter] to launch the juProt server..."
echo
echo "🚀 Starting the juProt application server..."
echo "   Waiting for server to be ready... (this may take a few seconds)"
echo "   When ready, access the application at http://localhost:8888"
echo "   Press [Ctrl+C] in this terminal to stop the server."
echo

# Run the container in INTERACTIVE mode (--rm -it)
# This shows logs directly and cleans up automatically on exit.
docker run --rm -it -p 8888:8080 --name $CONTAINER_NAME $IMAGE_NAME

echo
echo "✅ Server stopped successfully."
