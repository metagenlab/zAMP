#!/bin/bash

VERSION=$1
echo $1
CPU=$2
echo $2
GITHUBAT=$3
echo $3

docker build https://$GITHUBAT@github.com/metagenlab/microbiome16S_pipeline.git#$VERSION \
    -t metagenlab/amplicon_pipeline:$VERSION \
    -f ./Dockerfile \
    --build-arg GITHUB_AT=$GITHUBAT \
    --build-arg TEST_CPU=$CPU

docker build https://$GITHUBAT@github.com/metagenlab/microbiome16S_pipeline.git#$VERSION \
    -f ./Dockerfile.validation \
    --build-arg VERSION=$VERSION \
    --build-arg TEST_CPU=$CPU

docker push metagenlab/amplicon_pipeline:$VERSION
