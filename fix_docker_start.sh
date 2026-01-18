#!/bin/bash
set -e

# Define variables
NETWORK_NAME="crossbarv2-kg_default"
DATA_DIR="/GenSIvePFS/users/clzeng/docker_data"
WORKSPACE_DIR="/GenSIvePFS/users/clzeng/workspace/CROssBARv2-KG"

# Create network if it doesn't exist
docker network create $NETWORK_NAME 2>/dev/null || true

echo "Starting import container..."
# Run import container
# Note: mimicking the import service in docker-compose.yml
docker run --rm \
  --name import \
  --network $NETWORK_NAME \
  -e NEO4J_AUTH=none \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e FILL_DB_ON_STARTUP="yes" \
  -v "$DATA_DIR:/data" \
  -v "$WORKSPACE_DIR/scripts/import_local.sh:/scripts/import_local.sh" \
  -v "$WORKSPACE_DIR/biocypher-out:/biocypher-out" \
  -v "$WORKSPACE_DIR/biocypher-out:/GenSIvePFS/users/clzeng/workspace/CROssBARv2-KG/biocypher-out" \
  neo4j:5.26.19-enterprise \
  /bin/bash /scripts/import_local.sh

echo "Import container finished."

echo "Starting deploy container..."
# Run deploy container
# Note: mimicking the deploy service in docker-compose.yml
docker run -d \
  --name deploy \
  --network $NETWORK_NAME \
  -p 127.0.0.1:7474:7474 \
  -p 7687:7687 \
  -e NEO4J_dbms_security_auth__enabled="false" \
  -e NEO4J_dbms_databases_default__to__read__only="false" \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_PLUGINS='["apoc"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="apoc.*" \
  -v "$DATA_DIR:/data" \
  neo4j:5.26.19-enterprise

echo "Deploy container started."
