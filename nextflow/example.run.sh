#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-08-07 17:10:29
# @DESCRIPTION:

# Number of input parameters
param=$#
nextflow run rmatsTurbo.nf --config nextflow.config -resume
