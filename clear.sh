#!/bin/bash

# Function to remove a subfolder
remove_subfolder() {
    folder_name=$1
    if [ -d "$folder_name" ]; then
        rm -r "$folder_name"
        echo "Removed $folder_name/"
    else
        echo "Folder $folder_name does not exist."
    fi
}

# If the argument is "ALL", remove all subfolders
if [ "$1" == "ALL" ]; then
    remove_subfolder "output"
    remove_subfolder "figs"
else
    # Remove the specified subfolder from both output/ and figs/
    remove_subfolder "output/$1"
    remove_subfolder "figs/$1"
fi