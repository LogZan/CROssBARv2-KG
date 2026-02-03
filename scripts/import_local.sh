#!/bin/bash
# Import from local biocypher-out directory (skipping the build step)
# Uses the latest directory in /biocypher-out

sleep 2

# # Find the output directory (prefer the collective import call folder if it exists)
# if [ -d "/biocypher-out/crossbar_all_organisms_admin_import_call/" ]; then
#     LATEST_DIR="/biocypher-out/crossbar_all_organisms_admin_import_call/"
# else
#     LATEST_DIR=$(ls -td /biocypher-out/*/ 2>/dev/null | head -1)
# fi

# Use the output directory
LATEST_DIR="/biocypher-out/v1.2/"

if [ -z "$LATEST_DIR" ]; then
    echo "ERROR: No output directories found in /biocypher-out/"
    exit 1
fi

IMPORT_SCRIPT="${LATEST_DIR}neo4j-admin-import-call.sh"

echo "Using output directory: $LATEST_DIR"
echo "Looking for import script: $IMPORT_SCRIPT"

if [ -f "$IMPORT_SCRIPT" ]; then
    echo "Found import script, executing with bash..."
    # Execute the script with bash instead of making it executable
    bash "$IMPORT_SCRIPT"
else
    echo "ERROR: Import script not found: $IMPORT_SCRIPT"
    ls -la "$LATEST_DIR"
    exit 1
fi

neo4j start
sleep 10
neo4j stop
