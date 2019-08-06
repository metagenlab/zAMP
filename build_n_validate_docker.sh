#!/bin/bash

_VERSION=$0
_CPU=$1
_GITHUB_AT=$2

docker build -f ./Dockerfile . --build-arg GITHUB_AT=$_GITHUB_AT --build-arg VERSION=$_VERSION -t metagenlab/amplicon_pipeline:$_VERSION

docker build -f ./validation.Docker . --build-arg VERSION=$_VERSION --build-arg GITHUB_AT=$_GITHUB_AT TEST_CPU=$_CPU
