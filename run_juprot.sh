#!/bin/bash

# --- PART 1: Banner & Info (Restored "& Team") ---
echo "========================================================================"
echo "          juProt v1.0 - Comparative Protein-Ligand Analyzer"
echo
echo "           Developed by: Dr. Benedict Christopher Paul & Team"
echo "           Website: http://www.drpaul.cc"
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

# Check if the user has loaded the application image
if ! docker image inspect $IMAGE_NAME > /dev/null 2>&1; then
    echo "🔎 juProt image not found. Loading it from juprot.tar..."
    if [ -f "juprot.tar" ]; then
        docker load < juprot.tar
        echo "✅ Image loaded successfully."
    else
        echo "❌ Error: juprot.tar file not found. Make sure it's in the same folder."
        exit 1
    fi; echo
fi

# Check if the container is already running
if [ "$(docker ps -q -f name=$CONTAINER_NAME)" ]; then
    echo "🟢 juProt is already running. Access it at http://localhost:8888"
    exit 0
fi

# --- PART 3: LAUNCH ---
read -p "Press [Enter] to launch the juProt server..."
echo
echo "🚀 Starting the juProt application server..."
echo "   When ready, access the application at http://localhost:8888"
echo "   Press [Ctrl+C] in this terminal (once) to stop the server."
echo

# Use 'exec' to ensure Ctrl+C works on the first press.
# Update port mapping to 8888:8888 to match the Dockerfile.
exec docker run --rm -it -p 8888:8888 --name $CONTAINER_NAME $IMAGE_NAME