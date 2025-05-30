#!/bin/bash

# Number of screens to create
NUM_SCREENS=8

# Create screens
for i in $(seq 1 $NUM_SCREENS); do
    # Create a unique screen name
    SCREEN_NAME="s$i"
    
    # Create detached screen session
    screen -dmS "$SCREEN_NAME"
    
    echo "Created screen session: $SCREEN_NAME"
done

echo "Successfully created $NUM_SCREENS screen sessions"
echo "List of all screen sessions:"
screen -ls
