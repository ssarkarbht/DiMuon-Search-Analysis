#!/bin/python

'''
Author: Sourav Sarkar
Date: November 5, 2022
Email: ssarkar1@ualberta.ca
Description: This script processes all final MC datasets and produces
	pickle out file that are used in the final fit.
'''

import time
import json
from optparse import OptionParser
import sys
import os
cwd = os.getcwd()
moduledir = cwd+'/../modules/'
if moduledir not in sys.path:
    sys.path.append(moduledir)


from analysis_process import *

start_time = time.time()

parser=OptionParser()
parser.add_option("-c", "--ConfigFile", dest="CONFIG", type=str)
parser.add_option("-a", "--CreateAsimov", action = 'store_true', dest="ASIMOV", default=False)

(options, args) = parser.parse_args()

with open(options.CONFIG) as f:
	config = json.load(f)
asimov = options.ASIMOV

#******* Common Args *******
xsdir  = config["XSDir"]
nsqdir = config["NSQDir"]
charm_xs = config["CharmXS"]
band = ["UP", "DN"]

#Initialize placeholder list of analysis/process objects
sig_holder = []
bkg_holder = []
sigsys_holder = []
bkgsys_holder = []

nom_dict = config["Input"]["Nominal"]

#============== Process CCDIS single Muon ============================
infiles = [nom_dict["mc_loc"]+f for f in nom_dict["ccdis_files"]]
ascores = [nom_dict["score_loc"]+f for f in nom_dict["ccdis_Ascores"]]
bscores = [nom_dict["score_loc"]+f for f in nom_dict["ccdis_Bscores"]]
norm = config["Normalization"]["Nominal"]["ccdis"]

print ("---------- Process [1/9] : NuMu CCDIS single muon background")
# Process Nominal
nom_ccdis = process_dataset(infiles, ascores, bscores, 'ccdis', norm, xsdir, nsqdir)

# Process Weight based systematics
for w in config["WeightSystematics"]:
	for b in band:
		sysname = w+b
		print ("Processing weight systematics : ", sysname)
		sys = process_sys(infiles, norm, 'ccdis', sysname, xsdir, nsqdir)
		nom_ccdis.add_systematics_weight(sys[1], sysname)

bkg_holder.append(nom_ccdis)

#============== Process Charm dimuon =================================
infiles = [nom_dict["mc_loc"]+f for f in nom_dict["charm_files"]]
ascores = [nom_dict["score_loc"]+f for f in nom_dict["charm_Ascores"]]
bscores = [nom_dict["score_loc"]+f for f in nom_dict["charm_Bscores"]]
norm = config["Normalization"]["Nominal"]["charm"]

print ("---------- Process [2/9] : NuMu Charm DiMuon signal")
# Process Nominal
nom_charm = process_dataset(infiles, ascores, bscores, 'charmNom', norm, xsdir, nsqdir, xsdata=charm_xs)

# Process Weight based systematics
for w in config["WeightSystematics"]:
	for b in band:
		sysname = w+b
		print ("Processing weight systematics : ", sysname)
		sys = process_sys(infiles, norm, 'charm', sysname, xsdir, nsqdir, xsdata=charm_xs)
		nom_charm.add_systematics_weight(sys[1], sysname)

sig_holder.append(nom_charm)

#============== Process Trident dimuon ===============================
infiles = [nom_dict["mc_loc"]+f for f in nom_dict["trident_files"]]
ascores = [nom_dict["score_loc"]+f for f in nom_dict["trident_Ascores"]]
bscores = [nom_dict["score_loc"]+f for f in nom_dict["trident_Bscores"]]
norm = config["Normalization"]["Nominal"]["trident"]

print ("---------- Process [3/9] : NuMu Trident DiMuon signal")
nom_trident = process_dataset(infiles, ascores, bscores, 'trident', norm, xsdir, nsqdir)
sig_holder.append(nom_trident)

#============== Process NuTau Background =============================
bkg_dict = config["Input"]["AddBackground"]

infiles = [bkg_dict["mc_loc"] + bkg_dict["NuTau"]["file"]]
ascores = [bkg_dict["score_loc"] + bkg_dict["NuTau"]["Ascore"]]
bscores = [bkg_dict["score_loc"] + bkg_dict["NuTau"]["Bscore"]]
norm = config["Normalization"]["AddBackground"]["NuTau"]

print ("---------- Process [4/9] : NuTau CCDIS single muon background")
nutau_bkg = process_dataset(infiles, ascores, bscores, 'ccdis', norm, xsdir, nsqdir)
bkg_holder.append(nutau_bkg)

#============== Process NuE Background ===============================
infiles = [bkg_dict["mc_loc"] + bkg_dict["NuE"]["file"]]
ascores = [bkg_dict["score_loc"] + bkg_dict["NuE"]["Ascore"]]
bscores = [bkg_dict["score_loc"] + bkg_dict["NuE"]["Bscore"]]
norm = config["Normalization"]["AddBackground"]["NuE"]

print ("---------- Process [5/9] : NuE CCDIS single muon background")
nue_bkg = process_dataset(infiles, ascores, bscores, 'ccdis', norm, xsdir, nsqdir)
bkg_holder.append(nue_bkg)

#============== Process Corsika Background ===========================
infiles = [bkg_dict["mc_loc"] + bkg_dict["Corsika"]["file"]]
ascores = [bkg_dict["score_loc"] + bkg_dict["Corsika"]["Ascore"]]
bscores = [bkg_dict["score_loc"] + bkg_dict["Corsika"]["Bscore"]]
norm = config["Normalization"]["AddBackground"]["Corsika"]

print ("---------- Process [6/9] : CORSIKA Atmospheric muon background")
corsika = process_corsika(infiles, ascores, bscores, norm)
bkg_holder.append(corsika)

#============== Process DOM eff. MC sets =============================
print ("---------- Process [7/9] : DOM Eff. MC sets signal + background")
sys_dict = config["Input"]["Systematics"]
floc = config["Input"]["Systematics"]["mc_loc"]
sloc = config["Input"]["Systematics"]["score_loc"]

# Signal Dom Eff.
norm = config["Normalization"]["Systematics"]["charm"]["dom"]
infiles = [floc+f for f in sys_dict["charm_files"]["dom"]]
ascores = [sloc+f for f in sys_dict["charm_Ascores"]["dom"]]
bscores = [sloc+f for f in sys_dict["charm_Bscores"]["dom"]]

for i in range(len(infiles)):
	f = [infiles[i]]
	a = [ascores[i]]
	b = [bscores[i]]
	proc = process_dataset(f, a, b, 'charmSys', norm, xsdir, nsqdir, xsdata=charm_xs)
	sigsys_holder.append(proc)

#Background DOM eff
infiles = [floc+f for f in sys_dict["ccdis_files"]["dom"]]
ascores = [sloc+f for f in sys_dict["ccdis_Ascores"]["dom"]]
bscores = [sloc+f for f in sys_dict["ccdis_Bscores"]["dom"]]
norm = config["Normalization"]["Systematics"]["ccdis"]["dom"]

for i in range(len(infiles)):
	f = [infiles[i]]
	a = [ascores[i]]
	b = [bscores[i]]
	n = norm[i]
	proc = process_dataset(f, a, b, 'ccdis', n, xsdir, nsqdir)
	bkgsys_holder.append(proc)

#============== Process Hole Ice MC sets =============================
print ("---------- Process [8/9] : Hole Ice MC sets signal + background")
# Signal Hole Ice
norm = config["Normalization"]["Systematics"]["charm"]["hole"]
infiles = [floc+f for f in sys_dict["charm_files"]["hole"]]
ascores = [sloc+f for f in sys_dict["charm_Ascores"]["hole"]]
bscores = [sloc+f for f in sys_dict["charm_Bscores"]["hole"]]

for i in range(len(infiles)):
	f = [infiles[i]]
	a = [ascores[i]]
	b = [bscores[i]]
	proc = process_dataset(f, a, b, 'charmSys', norm, xsdir, nsqdir, xsdata=charm_xs)
	sigsys_holder.append(proc)

#Background Hole Ice
infiles = [floc+f for f in sys_dict["ccdis_files"]["hole"]]
ascores = [sloc+f for f in sys_dict["ccdis_Ascores"]["hole"]]
bscores = [sloc+f for f in sys_dict["ccdis_Bscores"]["hole"]]
norm = config["Normalization"]["Systematics"]["ccdis"]["hole"]

for i in range(len(infiles)):
	f = [infiles[i]]
	a = [ascores[i]]
	b = [bscores[i]]
	n = norm[i]
	proc = process_dataset(f, a, b, 'ccdis', n, xsdir, nsqdir)
	bkgsys_holder.append(proc)

#============== Process Bulk Ice MC sets =============================
print ("---------- Process [9/9] : Bulk Ice MC sets signal + background")
# Signal Bulk Ice
norm = config["Normalization"]["Systematics"]["charm"]["bulk"]
infiles = [floc+f for f in sys_dict["charm_files"]["bulk"]]
ascores = [sloc+f for f in sys_dict["charm_Ascores"]["bulk"]]
bscores = [sloc+f for f in sys_dict["charm_Bscores"]["bulk"]]

for i in range(len(infiles)):
	f = [infiles[i]]
	a = [ascores[i]]
	b = [bscores[i]]
	proc = process_dataset(f, a, b, 'charmSys', norm, xsdir, nsqdir, xsdata=charm_xs)
	sigsys_holder.append(proc)

# Background Bulk Ice
infiles = [floc + sys_dict["ccdis_files"]["bulk"][0]]
ascores = [sloc + sys_dict["ccdis_Ascores"]["bulk"][0]]
bscores = [sloc + sys_dict["ccdis_Bscores"]["bulk"][0]]
norm = config["Normalization"]["Systematics"]["ccdis"]["bulk"]

proc = process_dataset(infiles, ascores, bscores, 'ccdis', norm, xsdir, nsqdir)
bkgsys_holder.append(proc)



#=============== Get final event counts in  Cut regions
reglist = config["CutRegions"]
wtsyslist = config["WeightSystematics"]

#set up output directory
outdir = config["Output"]["outdir"]
if outdir is None:
	print (f"Output directory not specified in config. Will be storing output in current directory {cwd}")
	outdir = cwd+'/'
#******************** Burn Sample Expectation
livetime = config["Livetimes"]["burnsample"]
print ("Computing MC expectations for 0.76 Years (Burn Sample Livetime): ")

for key, value in reglist.items():
	print (f"---------------Event Expectation in {key}-----------------")
	tot, bkg, sig = get_MCCounts(sig_holder, bkg_holder,
			 sigsys_holder, bkgsys_holder, wtsyslist, value, livetime)
	print (f"Total (S+B)    : {tot[0]}")
	print (f"Signal (S)     : {sig[0]}")
	print (f"Background (B) : {bkg[0]}")

outfile = outdir+config["Output"]["outfiles"]["mc_burn"]
print (f"Storing expectations in pickle file: {outfile}")

generate_mcPickle(reglist, sig_holder, bkg_holder, sigsys_holder, bkgsys_holder, wtsyslist,
		livetime, outfile, asimovfile=None)

#******************** 10.67 Years Expectation
livetime = config["Livetimes"]["fullsample"]
print ("Computing MC expectations for 10.67 Years (Full Sample Livetime): ")
for key, value in reglist.items():
	print (f"---------------Event Expectation in {key}-----------------")
	tot, bkg, sig = get_MCCounts(sig_holder, bkg_holder,
			 sigsys_holder, bkgsys_holder, wtsyslist, value, livetime)
	print (f"Total (S+B)    : {tot[0]}")
	print (f"Signal (S)     : {sig[0]}")
	print (f"Background (B) : {bkg[0]}")

outfile = outdir+config["Output"]["outfiles"]["mc_full"]
print (f"Storing expectations in pickle file: {outfile}")

outfile_asimov = None
if asimov:
	outfile_asimov = outdir+config["Output"]["outfiles"]["mc_asimov"]
	print (f"Will be generating asimov dataset in {outfile_asimov}")

generate_mcPickle(reglist, sig_holder, bkg_holder, sigsys_holder, bkgsys_holder, wtsyslist,
                livetime, outfile, asimovfile=outfile_asimov)

#-------------------------- End of processing ---------------------------------------------
end_time = time.time()
duration = (end_time-start_time)/60.
print (f"Took {duration} minutes to run! GoodBye!")
