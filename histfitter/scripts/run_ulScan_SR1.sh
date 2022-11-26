#!/bin/bash

seed=100
#go to histfitter install directory
cd /usr/local/HistFitter/
#set up histfitter
source setup.sh
wait
#go to histfitter analysis directory
cd /DiMuonAnalysis/histfitter

#Perform the parameter fit
HistFitter.py -t -w -f scripts/ulscan_sr1_asimov.py
wait

#Perform Upper Limit scan for signal strength parameter
HistFitter.py -F excl -l -s $seed scripts/ulscan_sr1_asimov.py
