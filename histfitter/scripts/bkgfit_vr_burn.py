"""
Author: Sourav Sarkar
Email: ssarkar1@ualberta.ca
Date : November 08, 2022
Description:
	This Histfitter configuration file performs a
	background-only fit considering CR and VRs
	with only background event counts as the input
	from MC PDFs (+ uncertainity) and observed events
	in CR, VRs
"""

#import stuff
from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt
import pickle
import os

from ROOT import gROOT
#gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
#ROOT.SetAtlasStyle()

#### Define data inputs from pickle file
datadir = "/DiMuonAnalysis/process_data/"

infile_mc = datadir + "MCBurnInput.pkl"
with open(infile_mc, 'rb') as f:
	inputs = pickle.load(f)

infile_data = datadir + "DataBurn.pkl"
with open(infile_data, 'rb') as f:
	data = pickle.load(f)

#------> Inputs from Background MC template
iCR  = inputs['CR']["Background"]
iVR1 = inputs['VR1']["Background"]
iVR2 = inputs['VR2']["Background"]
iVR3 = inputs['VR3']["Background"]
iVR4 = inputs['VR4']["Background"]

#%%%%%%%%%%%%%%%%%%%%%%% Analysis Configuration %%%%%%%%%%%%%%%%%%%%%%%%
if not os.path.exists("results"):
	os.makedirs("results")

# Give the analysis a name
configMgr.analysisName = "BKGFitBurn"
configMgr.outputFileName = f"results/{configMgr.analysisName}s_Output.root"

# Setting the parameters of the hypothesis test
configMgr.doExclusion=True # True=exclusion, False=discovery
configMgr.nTOYs=10000
configMgr.calculatorType=0 # 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=3   # 3=one-sided profile likelihood test statistic (LHC default)
configMgr.nPoints=20       # number of values scanned of signal-strength for upper-limit determination of signal strength.
configMgr.writeXML = False

#%%%%%%%%%%%%%%%%%%%%%%% Define Channels, Samples, Systematics %%%%%%%%%
# Define cut Regions/Channels
configMgr.cutsDict["CR"] = "1." # "1." denotes the cut is manually hadnled
configMgr.cutsDict["VR1"] = "1."
configMgr.cutsDict["VR2"] = "1."
configMgr.cutsDict["VR3"] = "1."
configMgr.cutsDict["VR4"] = "1."

# Define weights
configMgr.weights = "1."

# Define MC Samples
bkgSample = Sample("Bkg", kGreen-9)
bkgSample.setNormFactor("mu_BG",1.,0.,5.)
bkgSample.setStatConfig(True)

# Add Nominal counts to CR, VRs
# "cuts" for cut&count, 0.5 for lower bin edge (dummy for 1 bin)
bkgSample.buildHisto([iCR["NEvents"]],"CR","cuts",0.5)
bkgSample.buildHisto([iVR1["NEvents"]],"VR1","cuts",0.5) 
bkgSample.buildHisto([iVR2["NEvents"]],"VR2","cuts",0.5) 
bkgSample.buildHisto([iVR3["NEvents"]],"VR3","cuts",0.5) 
bkgSample.buildHisto([iVR4["NEvents"]],"VR4","cuts",0.5)

# Add MC Stat. Errors to CR, VRs
bkgSample.buildStatErrors([iCR["MCUnc"]],"CR","cuts")
bkgSample.buildStatErrors([iVR1["MCUnc"]],"VR1","cuts")
bkgSample.buildStatErrors([iVR2["MCUnc"]],"VR2","cuts")
bkgSample.buildStatErrors([iVR3["MCUnc"]],"VR3","cuts")
bkgSample.buildStatErrors([iVR4["MCUnc"]],"VR4","cuts")

# Define data samples
dataSample = Sample("Data", kBlack)
dataSample.setData()

dataSample.buildHisto([data["CR"]],"CR","cuts",0.5)
dataSample.buildHisto([data["VR1"]],"VR1","cuts",0.5)
dataSample.buildHisto([data["VR2"]],"VR2","cuts",0.5)
dataSample.buildHisto([data["VR3"]],"VR3","cuts",0.5)
dataSample.buildHisto([data["VR4"]],"VR4","cuts",0.5)


# Define MC Sys. Errors to CR, VRs
sysmethod = "userNormHistoSys"

# correlated systematic between backgrounds in different regions 1 +- relative uncertainties
sys1CR  = Systematic("corrsys",configMgr.weights, [iCR["CorrUp"]],[iCR["CorrDn"]], "user", sysmethod)
sys1VR1 = Systematic("corrsys",configMgr.weights, [iVR1["CorrUp"]],[iVR1["CorrDn"]], "user", sysmethod)
sys1VR2 = Systematic("corrsys",configMgr.weights, [iVR2["CorrUp"]],[iVR2["CorrDn"]], "user", sysmethod)
sys1VR3 = Systematic("corrsys",configMgr.weights, [iVR3["CorrUp"]],[iVR3["CorrDn"]], "user", sysmethod)
sys1VR4 = Systematic("corrsys",configMgr.weights, [iVR4["CorrUp"]],[iVR4["CorrDn"]], "user", sysmethod)

# uncorrelated systematic in different regions
sys2CR  = Systematic("uncorrCR",configMgr.weights, 1.+iCR["UnCorr"], 1.-iCR["UnCorr"], "user", "userOverallSys")
sys2VR1 = Systematic("uncorrVR1",configMgr.weights, 1.+iVR1["UnCorr"], 1.-iVR1["UnCorr"], "user", "userOverallSys")
sys2VR2 = Systematic("uncorrVR2",configMgr.weights, 1.+iVR2["UnCorr"], 1.-iVR2["UnCorr"], "user", "userOverallSys")
sys2VR3 = Systematic("uncorrVR3",configMgr.weights, 1.+iVR3["UnCorr"], 1.-iVR3["UnCorr"], "user", "userOverallSys")
sys2VR4 = Systematic("uncorrVR4",configMgr.weights, 1.+iVR4["UnCorr"], 1.-iVR4["UnCorr"], "user", "userOverallSys")

#%%%%%%%%%%%%%% Build the model: fitConfig, add channels, add samples, add systematics

# Define fitConfig
bkg_fit = configMgr.addFitConfig("BkgFit")

# Define measurement
meas = bkg_fit.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=0.0)
meas.addPOI("mu_Sig")

meas.addParamSetting("Lumi",True,1)

# Add the channels
RegC  = bkg_fit.addChannel("cuts", ["CR"], 1,0.5,1.5)
RegV1 = bkg_fit.addChannel("cuts", ["VR1"], 1,0.5,1.5)
RegV2 = bkg_fit.addChannel("cuts", ["VR2"], 1,0.5,1.5)
RegV3 = bkg_fit.addChannel("cuts", ["VR3"], 1,0.5,1.5)
RegV4 = bkg_fit.addChannel("cuts", ["VR4"], 1,0.5,1.5)

# Add the samples
bkg_fit.addSamples([bkgSample, dataSample])

# Add the systematics
RegC.getSample("Bkg").addSystematic(sys1CR)
RegC.getSample("Bkg").addSystematic(sys2CR)

RegV1.getSample("Bkg").addSystematic(sys1VR1)
RegV1.getSample("Bkg").addSystematic(sys2VR1)

RegV2.getSample("Bkg").addSystematic(sys1VR2)
RegV2.getSample("Bkg").addSystematic(sys2VR2)

RegV3.getSample("Bkg").addSystematic(sys1VR3)
RegV3.getSample("Bkg").addSystematic(sys2VR3)

RegV4.getSample("Bkg").addSystematic(sys1VR4)
RegV4.getSample("Bkg").addSystematic(sys2VR4)


# Define the channel labels
bkg_fit.addBkgConstrainChannels([RegC])
bkg_fit.addValidationChannels([RegV1, RegV2, RegV3, RegV4])


# These lines are needed for the user analysis to run
# Make sure file is re-made when executing HistFactory
if configMgr.executeHistFactory:
	if os.path.isfile(f"data/{configMgr.analysisName}.root"):
		os.remove(f"data/{configMgr.analysisName}.root")


