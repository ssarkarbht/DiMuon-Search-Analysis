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
HistFitter.py -t -w -f scripts/discfit_sr2_asimov.py
wait

#get the systematics table
SysTable.py -o results/disc_SR2Syst.tex -c CR,SR -s Bkg,Sig -w results/SR2_DiscFit/SPlusB_combined_NormalMeasurement_model_afterFit.root -%
wait

#get the fitted event counts in each region (with error)
YieldsTable.py -o results/disc_SR2Yield.tex -c CR,SR -s Bkg,Sig -w results/SR2_DiscFit/SPlusB_combined_NormalMeasurement_model_afterFit.root -C "SR2 Fit on Asimov Dataset (S+B)"
wait

#Perform Discovery fit hypothesis test
HistFitter.py -F disc -z -s $seed scripts/discfit_sr2_asimov.py
