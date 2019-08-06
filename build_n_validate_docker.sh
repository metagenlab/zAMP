#!/bin/bash

VERSION=$1
CPU=$2
GITHUBAT=$3

docker build https://$GITHUBAT@github.com/metagenlab/microbiome16S_pipeline.git#$VERSION -t metagenlab/amplicon_pipeline:$VERSION -f ./Dockerfile . --build-arg GITHUB_AT=$_GITHUB_AT --build-arg VERSION=$_VERSION

docker build -f ./Dockerfile . --build-arg GITHUB_AT=$_GITHUB_AT --build-arg VERSION=$_VERSION -t $_VERSION

docker build -f ./validation.Docker . --build-arg VERSION=$_VERSION --build-arg GITHUB_AT=$_GITHUB_AT TEST_CPU=$_CPU
