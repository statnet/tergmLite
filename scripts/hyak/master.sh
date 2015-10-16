#!/bin/bash

qsub -q batch runsim.tlite.sh
qsub -q batch runsim.tbig.sh
