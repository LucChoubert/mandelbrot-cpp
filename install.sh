#!/bin/bash

USER=`whoami`
TARGET_DIR="/home/$USER/.local/bin"
FILE="Mandelbrot"

if [[ -d "$TARGET_DIR" ]]
then
    if [[ -f "$FILE" ]]
    then
        echo "Installing $FILE on your local bin directory $TARGET_DIR"
        cp $FILE $TARGET_DIR
        echo "Done!"
    else
        echo "$FILE does not exist. Exiting..."
        exit 1
    fi
else
    echo "No local bin directory $TARGET_DIR. Exiting..."
    exit 1
fi


