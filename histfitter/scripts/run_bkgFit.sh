#!/bin/bash

seed=100
#go to histfitter install directory
cd /usr/local/HistFitter/
#set up histfitter
source setup.sh
wait
#go to histfitter analysis directory
cd /DiMuonAnalysis/histfitter

echo `pwd`
#Perform the parameter fit
HistFitter.py -t -w -f -F bkg -s $seed scripts/bkgfit_vr_burn.py
wait

#get the systematics table
SysTable.py -o results/burnSyst.tex -c CR,VR1,VR2,VR3,VR4 -s Bkg -w results/BKGFitBurn/BkgFit_combined_NormalMeasurement_model_afterFit.root -%
wait

#get the fitted event counts in each region (with error)
YieldsTable.py -o results/burnYield.tex -c CR,VR1,VR2,VR3,VR4 -s Bkg -w results/BKGFitBurn/BkgFit_combined_NormalMeasurement_model_afterFit.root -C "Background-only fit on Burn Sample"


