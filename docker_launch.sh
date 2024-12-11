#!/bin/bash
TOOL_NAME=$(Rscript -e "cat(desc::desc_get_field('Package'))")
TOOL_VERSION=$(Rscript -e "cat(desc::desc_get_field('Version'))")
DOCKER_NAME=$(Rscript -e "cat(desc::desc_get_field('Package'))")
DOCKER_VERSION="1.0.0"

docker build -t ${DOCKER_NAME}:${DOCKER_VERSION} --build-arg TOOL_NAME=${TOOL_NAME} --build-arg TOOL_VERSION=${TOOL_VERSION} dev
