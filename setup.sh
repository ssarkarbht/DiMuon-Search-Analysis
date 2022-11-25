#!/bin/bash

#Setup script to run data processing scripts

#Load the cvmfs environment
echo "Loading CVMFS environment ..."
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

#Add the library path needed for running nuSQuIDs
echo "Adding nuSQuIDs library path ..."
export LD_LIBRARY_PATH=/data/user/ssarkar/I3/LW/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/data/user/ssarkar/I3/LW/lib:$PYTHONPATH

#Add current directory as the parent analysis directory
echo "Setting Analysis driectory variable name ..."
export ANALYSIS_DIR=$(echo `pwd`)

echo "Done."

