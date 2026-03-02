#!/bin/bash

CONTAINER_NAME="juprot-container"

if [ ! "$(docker ps -q -f name=$CONTAINER_NAME)" ]; then
    echo "juProt is not currently running."
    # Clean up any old, stopped containers with that name
    if [ "$(docker ps -aq -f status=exited -f name=$CONTAINER_NAME)" ]; then
        docker rm $CONTAINER_NAME > /dev/null
    fi
    exit 0
fi

echo "Found juProt running in the background. Stopping it now..."
docker stop $CONTAINER_NAME > /dev/null
docker rm $CONTAINER_NAME > /dev/null
echo "✅ juProt has been stopped."
