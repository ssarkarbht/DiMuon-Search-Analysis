"""
Author: Sourav Sarkar
Email: ssarkar1@ualberta.ca
Date : November 08, 2022
Description:
        This Histfitter configuration file performs a
        signal fit considering CR and SR
        with event counts as the input
        from MC PDFs (+ uncertainities) and
	fake/asimov as observed data
"""
#import stuff
from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt
import pickle

from ROOT import gROOT
#gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
#ROOT.SetAtlasStyle()

#### Define data inputs from pickle file
datadir = "/DiMuonAnalysis/process_data/"

infile_mc = datadir+"MCFullInput.pkl"

with open(infile_mc, 'rb') as f:
        inputs = pickle.load(f)

infile_data = datadir+"MCAsimovInput.pkl"

with open(infile_data, 'rb') as f:
        asimov_data = pickle.load(f)

signalR = "SR1"
#******-------------- Background Sample Inputs ----------------
bkgCR = inputs['CR']["Background"]
bkgSR = inputs[signalR]["Background"]
#get ROI for calculating systematics in SRs
bkgROI = inputs['ROI']["Background"]

#******-------------- Signal Sample Inputs -------------------
sigCR = inputs['CR']["Signal"]
sigSR = inputs[signalR]["Signal"]
#get ROI for calculating systematics in SRs
sigROI = inputs['ROI']["Signal"]

#******-------------- Data Sample Inputs ---------------------
data_cr = asimov_data["CR"]["asimov_s+b"]
data_sr = asimov_data[signalR]["asimov_s+b"]

##########################
anaName = signalR+"_DiscFit"

# Setting the parameters of the hypothesis test
configMgr.doExclusion=False # True=exclusion, False=discovery
configMgr.nTOYs=5000
configMgr.calculatorType=0 # 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=2   # 3=one-sided profile likelihood test statistic (LHC default)

configMgr.writeXML = False

# Keep SRs also in background fit confuguration
configMgr.keepSignalRegionType = True
##########################

# Give the analysis a name
configMgr.analysisName = anaName
configMgr.outputFileName = f"results/{configMgr.analysisName}s_Output.root"

#%%%%%%%%%%%%%%%%%%%%%%% Define Channels, Samples, Systematics %%%%%%%%%
setStat = True

configMgr.cutsDict["CR"] = "1." # "1." denotes the cut is manually hadnled
configMgr.cutsDict["SR"] = "1."

# Define weights
configMgr.weights = "1."

# Define MC Samples
bkgSample = Sample("Bkg", kGreen-9)
bkgSample.setNormFactor("mu_BG",1.,0.,5.)
bkgSample.setStatConfig(setStat)

bkgSample.buildHisto([bkgCR["NEvents"]],"CR","cuts",0.5)
bkgSample.buildStatErrors([bkgCR["MCUnc"]],"CR", "cuts")

bkgSample.buildHisto([bkgSR["NEvents"]],"SR","cuts",0.5)
bkgSample.buildStatErrors([bkgSR["MCUnc"]],"SR", "cuts")

sigSample = Sample("Sig", kPink)
sigSample.setNormFactor("mu_Sig", 1.,0.,10.)
sigSample.setStatConfig(setStat)
sigSample.setNormByTheory()

sigSample.buildHisto([sigCR["NEvents"]], "CR", "cuts",0.5)
sigSample.buildStatErrors([sigCR["MCUnc"]], "CR", "cuts")

sigSample.buildHisto([sigSR["NEvents"]], "SR", "cuts",0.5)
sigSample.buildStatErrors([sigSR["MCUnc"]], "SR", "cuts")


dataSample = Sample("Data", kBlack)
dataSample.setData()
dataSample.buildHisto([data_cr],"CR","cuts",0.5)
dataSample.buildHisto([data_sr],"SR","cuts",0.5)

# Define Systematics
sysmethod = "userNormHistoSys"
sys1CRbkg  = Systematic("corrsys",configMgr.weights, [bkgCR["CorrUp"]],[bkgCR['CorrDn']], "user", sysmethod)
sys1CRsig  = Systematic("corrsys",configMgr.weights, [sigCR['CorrUp']],[sigCR['CorrDn']], "user", sysmethod)

sys1SRbkg  = Systematic("corrsys",configMgr.weights, [bkgROI['CorrUp']],[bkgROI['CorrDn']], "user", sysmethod)
sys1SRsig  = Systematic("corrsys",configMgr.weights, [sigROI['CorrUp']],[sigROI['CorrDn']], "user", sysmethod)

sys2CRbkg = Systematic("uncsysCRbkg",configMgr.weights, 1.+bkgCR['UnCorr'], 1.-bkgCR['UnCorr'], "user", "userOverallSys")
sys2CRsig = Systematic("uncsysCRsig",configMgr.weights, 1.+sigCR['UnCorr'], 1.-sigCR['UnCorr'], "user", "userOverallSys")

sys2SRbkg = Systematic("uncsysSRbkg",configMgr.weights, 1.+bkgROI['UnCorr'], 1.-bkgROI['UnCorr'], "user", "userOverallSys")
sys2SRsig = Systematic("uncsysSRsig",configMgr.weights, 1.+sigROI['UnCorr'], 1.-sigROI['UnCorr'], "user", "userOverallSys")

# Define top-level
ana = configMgr.addFitConfig("SPlusB")
ana.addSamples([bkgSample,sigSample,dataSample])
ana.setSignalSample(sigSample)

# Define measurement
meas = ana.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=1e-10)
meas.addPOI("mu_Sig")
meas.addParamSetting("Lumi",True,1)

# Add the channels
cReg = ana.addChannel("cuts",["CR"],1,0.5,1.5)
sReg = ana.addChannel("cuts",["SR"],1,0.5,1.5)

ana.addBkgConstrainChannels([cReg])
ana.addSignalChannels([sReg])

# Add the systematics
cReg.getSample("Bkg").addSystematic(sys1CRbkg)
cReg.getSample("Bkg").addSystematic(sys2CRbkg)
cReg.getSample("Sig").addSystematic(sys1CRsig)
cReg.getSample("Sig").addSystematic(sys2CRsig)

sReg.getSample("Bkg").addSystematic(sys1SRbkg)
sReg.getSample("Bkg").addSystematic(sys2SRbkg)
sReg.getSample("Sig").addSystematic(sys1SRsig)
sReg.getSample("Sig").addSystematic(sys2SRsig)

# These lines are needed for the user analysis to run
# Make sure file is re-made when executing HistFactory
if configMgr.executeHistFactory:
    if os.path.isfile(f"data/{configMgr.analysisName}.root"):
        os.remove(f"data/{configMgr.analysisName}.root")


