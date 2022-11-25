#!/bin/python

'''
Author: Sourav Sarkar
Date: November 2, 2022
Email: ssarkar1@ualberta.ca
Description: This script loads the analysis module
	and processes all the final MC/data samples
	to get arrays for plotting, extracting event
	counts, compute systematic variations
'''
#import modues
from glob import glob
import pickle
import sys
import os
from os.path import expandvars
cwd = os.getcwd()

moduledir = cwd+'/'
if moduledir not in sys.path:
    sys.path.append(moduledir)

from analysis_module import *

#wrapper function for processing all physics MC sets
def process_dataset(dfiles, safiles, sbfiles, dtype, norm,
                    xsdir, nsqdir, xsdata=None):
    #Initialize true energy and weight array
    enarr = np.array([])
    wtarr = np.array([])
    
    print (f"Processing {len(dfiles)} files... sit tight!")
    # if processing CCDIS Backgrounds
    if dtype=='ccdis':
        for run,f in enumerate(dfiles):
            print (f"Processing {f}")
            weight = Reweighting(f, xsdir, nsqdir)
    
            print ("\rUpdating cross section...(Charm Correction)", end="\r")
            weight.update_cross_section()
            print ("\rCalculating Total Flux...", end="\r")
            weight.get_total_flux()
    
            enarr = np.append(enarr, weight.nuen)
            wtarr = np.append(wtarr, weight.onew)
        wtarr /= norm
    # if processing nominal Charm signal MC
    elif dtype=='charmNom':
        for run,f in enumerate(dfiles):
            print (f"Processing {f}")
            weight = Reweighting(f, xsdir, nsqdir)
            runnum='Run'+'{:02d}'.format(run+3)
            xs = weight.load_xsdata(xsdata, runnum)
            #run 3,4,5 MC sets need the interaction weight correction
            if run<3:
                print ("\rUpdating Interaction weight...", end="\r")
                weight.update_interaction_weight(xs)

            print ("\rMaking charm xs pdfset correction (CT14nlo to CT18Anlo)...", end="\r")
            weight.charm_pdf_correction(xs)
            print ("\rDecoupling Propagation Weight...", end="\r")
            weight.undo_propagation_weight()
            print ("\rApplying total flux", end="\r")
            weight.get_total_flux()
    
            enarr = np.append(enarr, weight.nuen)
            wtarr = np.append(wtarr, weight.onew)
            
        wtarr /= norm
        
    elif dtype=='charmSys' or dtype=='taucharm':
        for run,f in enumerate(dfiles):
            print (f"Processing {f}")
            weight = Reweighting(f, xsdir, nsqdir)
            
            print ("\rMaking charm xs pdfset correction (CT14nlo to CT18Anlo)...", end="\r")
            corr_factor = weight.cross_section.get_pdfset_correction(weight.nuen, weight.nupdg)
            weight.onew *= corr_factor
            print ("\rDecoupling Propagation Weight...", end="\r")
            weight.undo_propagation_weight()
            print ("\rApplying total flux", end="\r")
            weight.get_total_flux()
    
            enarr = np.append(enarr, weight.nuen)
            wtarr = np.append(wtarr, weight.onew)
            
        wtarr /= norm

    
    elif dtype=='trident':
        #There two Trident dataset with different injected nuen energy range
        # To merge both datasets, need to take average of the two dataset weights
        # in the neutrino energy range (1TeV - 1PeV)
        EMIN = 1e3
        EMAX = 1e6
        for run,f in enumerate(dfiles):
            print (f"Processing {f}")
            weight = Reweighting(f, xsdir, nsqdir)
            print ("\rDecoupling Propagation Weight...", end="\r")
            weight.undo_propagation_weight()
            print ("\rApplying total flux", end="\r")
            weight.get_total_flux()
            
            tempen = weight.nuen
            tempwt = weight.onew
            
            avg_idx = np.where((tempen>=EMIN)&(tempen<=1e6))[0]
            
            tempwt[avg_idx] = tempwt[avg_idx]/2.0
            
            enarr = np.append(enarr, tempen)
            wtarr = np.append(wtarr, tempwt)
            
        wtarr /= norm
    #Now build the final array
    final_arr = AnalysisCutV2(dfiles, safiles, sbfiles, wtarr)
    return final_arr

#wrapper function for computing weight based systematics
def process_sys(dfiles, norm, dtype, systype,
                xsdir, nsqdir, xsdata=None):
    fluxunc = ['ConvNormUP', 'ConvNormDN',
              'ConvGammaUP', 'ConvGammaDN',
              'AstroNormUP', 'AstroNormDN',
              'AstroGammaUP', 'AstroGammaDN']
    xsunc = ['XSUP', 'XSDN']
    
    enarr = np.array([])
    wtarr = np.array([])
    
    print (f"Processing {len(dfiles)} files...")
    for run,f in enumerate(dfiles):
        print (f"Processing {f}")
        weight = Reweighting(f, xsdir, nsqdir)
        if dtype=='ccdis' and systype=='XSUP':
            print ("\rUpdating cross section... (Charm Correction) with +10% unc.", end="\r")
            weight.update_cross_section(sys='up')
            
        if dtype=='ccdis' and systype=='XSDN':
            print ("\rUpdating cross section... (Charm Correction) with -10% unc.", end="\r")
            weight.update_cross_section(sys='dn')
        
        if dtype=='ccdis' and systype not in xsunc:
            print ("\rUpdating cross section... (Charm Correction)", end="\r")
            weight.update_cross_section()

        if dtype=='charm':
            
            runnum='Run'+'{:02d}'.format(run+3)
            xs = weight.load_xsdata(xsdata, runnum)
            if run<3:
                print ("\rUpdating Charm Interaction weight correction...", end="\r")
                weight.update_interaction_weight(xs)

            print ("\rMaking charm xs pdfset correction (CT14nlo to CT18Anlo)...", end="\r")
            weight.charm_pdf_correction(xs)

            if  systype=='XSUP':
                print ("\rUpdating Charm cross section... with +10% unc.", end="\r")
                weight.charm_xs_sys(xs, 'up')
                
            if systype=='XSDN':
                print ("\rUpdating Charm cross section... with -10% unc.", end="\r")
                weight.charm_xs_sys(xs, 'dn')
                
            print ("\rDecoupling Propagation Weight...", end="\r")
            weight.undo_propagation_weight()

        if systype=='ConvNormUP':
            print ("\rCalculating Total Flux... with +30% Atm. Conv. Normalization", end="\r")
            weight.get_total_flux(conv_norm='up')
            
        if systype=='ConvNormDN':
            print ("\rCalculating Total Flux... with -30% Atm. Conv. Normalization", end="\r")
            weight.get_total_flux(conv_norm='dn')
            
        if systype=='ConvGammaUP':
            print ("\rCalculating Total Flux... with +0.01 Atm. Conv. Gamma shift", end="\r")
            weight.get_total_flux(conv_gamma='up')
            
        if systype=='ConvGammaDN':
            print ("\rCalculating Total Flux... with -0.01 Atm. Conv. Gamma shift", end="\r")
            weight.get_total_flux(conv_gamma='dn')

        if systype=='AstroNormUP':
            print ("\rCalculating Total Flux... with +18% Astro Normalization", end="\r")
            weight.get_total_flux(astro_norm='up')
            
        if systype=='AstroNormDN':
            print ("\rCalculating Total Flux... with -18% Astro Normalization", end="\r")
            weight.get_total_flux(astro_norm='dn')
            
        if systype=='AstroGammaUP':
            print ("\rCalculating Total Flux... with +0.09 Astro Gamma shift", end="\r")
            weight.get_total_flux(astro_gamma='up')
            
        if systype=='AstroGammaDN':
            print ("\rCalculating Total Flux... with -0.09 Astro Gamma shift", end="\r")
            weight.get_total_flux(astro_gamma='dn')

        if systype not in fluxunc:
            print ("\rCalculating Total Flux...", end="\r")
            weight.get_total_flux()
    
        enarr = np.append(enarr, weight.nuen)
        wtarr = np.append(wtarr, weight.onew)

    wtarr /= norm
    return enarr, wtarr

#process the CORSIKA files separately as their weight definition is different
def process_corsika(dfiles, safiles, sbfiles, norm):
    wtarr = np.array([])
    for i,f in enumerate(dfiles):
        hf = h5.File(f,'r')
        wf = np.array(hf['WeightFactor'])/norm

        wtarr = np.append(wtarr, wf)

    final_arr = AnalysisCutV2(dfiles, safiles, sbfiles, wtarr)
    return final_arr

#wrapper function to process all signal+background MC (both nominal+systematics)
#sets and give the following values (inputs for the final fits):
# 1. Nominal event counts in the given cut region
# 2. MC Statistical Uncertainity (abs.) in the given cut region
# 3. Total uncorrelated systematic unc. (rel.) in the given region
# 4. Total correlated systematics unc. (rel.) lower variation
# 5. Total correlated systematics unc. (rel.) upper variation

def get_MCCounts(sigproc, bkgproc, sigsys, bkgsys, wtsys, cutreg, time):
    '''Args:
    sigproc (list): List of analysis objects containing the signal processes
    bkgproc (list): List of analysis objects containing the background processes
    sigsys (list): List of signal detector systematics analysis objects
    bkgsys (list): List of background detector systematics analysis objects
    wtsys (list): List of strings indicating the weight based systematics sets
    cutreg (tuple): cut region where we want to compute the outputs
    time (float): Total livetime for which we want to compute the outputs
    '''
    #apply analysis cuts
    for proc in sigproc+bkgproc+sigsys+bkgsys:
        proc.reset_cut()
        proc.apply_energyCut()
        proc.apply_zenithCut()
        if cutreg[0]=='boxcut':
            proc.apply_boxCut(cutreg[1][0], cutreg[1][1],
                             amax=cutreg[1][2],
                             bmax=cutreg[1][3])
        elif cutreg[0] in ['ellipse','hyperbola']:
            proc.apply_curveCut(cutreg[1][0], cutreg[1][1],
                               cutreg[0])

    #compute signals
    nsig = []
    nsig_mcerr = []
    TOTSIG = 0
    for proc in sigproc:
        nval = np.sum(proc.event_arr[proc.event_idx]['weight']*time)
        nerr = np.sqrt(np.sum((proc.event_arr[proc.event_idx]['weight']*time)**2))
        
        TOTSIG += nval
        nsig.append(nval)
        nsig_mcerr.append(nerr)
        
    #compute backgrounds
    nbkg = []
    nbkg_mcerr = []
    TOTBG = 0
    for proc in bkgproc:
        nval = np.sum(proc.event_arr[proc.event_idx]['weight']*time)
        nerr = np.sqrt(np.sum((proc.event_arr[proc.event_idx]['weight']*time)**2))
        
        TOTBG += nval
        nbkg.append(nval)
        nbkg_mcerr.append(nerr)
    
    #compute Signal + Background Det. Sys
    corrsigsysup = []
    corrsigsysdn = []
    corrbkgsysup = []
    corrbkgsysdn = []
    uncorrsigsys = []
    uncorrbkgsys = []
    for i in range(3):
        if i==2:#deal with 4 bulk variants separately
            #Signal
            n1 = np.sum(sigsys[4].event_arr[sigsys[4].event_idx]['weight']*time)
            n2 = np.sum(sigsys[5].event_arr[sigsys[5].event_idx]['weight']*time)
            n3 = np.sum(sigsys[6].event_arr[sigsys[6].event_idx]['weight']*time)
            n4 = np.sum(sigsys[7].event_arr[sigsys[7].event_idx]['weight']*time)
            nval = max(abs(nsig[0]-n1),
                       abs(nsig[0]-n2),
                       abs(nsig[0]-n3),
                       abs(nsig[0]-n4))
            uncorrsigsys.append(nval)
            
            #Background
            n = np.sum(bkgsys[4].event_arr[bkgsys[4].event_idx]['weight']*time)
            uncorrbkgsys.append(abs(nbkg[0]-n))
            
        elif i==1:#Hole ice systematics
            #signal
            n1 = np.sum(sigsys[2*i].event_arr[sigsys[2*i].event_idx]['weight']*time)
            n2 = np.sum(sigsys[2*i+1].event_arr[sigsys[2*i+1].event_idx]['weight']*time)
            nval = max(abs(nsig[0]-n1),
                       abs(nsig[0]-n2))
            uncorrsigsys.append(nval)

            #background
            n1 = np.sum(bkgsys[2*i].event_arr[bkgsys[2*i].event_idx]['weight']*time)
            n2 = np.sum(bkgsys[2*i+1].event_arr[bkgsys[2*i+1].event_idx]['weight']*time)
            nval = max(abs(nbkg[0]-n1),
                       abs(nbkg[0]-n2))
            uncorrbkgsys.append(nval)
        elif i==0:#DOM eff. systematics --> goes to correlated systematics
            #signal
            n1 = np.sum(sigsys[2*i].event_arr[sigsys[2*i].event_idx]['weight']*time)
            n2 = np.sum(sigsys[2*i+1].event_arr[sigsys[2*i+1].event_idx]['weight']*time)
            corrsigsysup.append(max(n1,n2)-nsig[0])
            corrsigsysdn.append(min(n1,n2)-nsig[0])
            
            #background
            n1 = np.sum(bkgsys[2*i].event_arr[bkgsys[2*i].event_idx]['weight']*time)
            n2 = np.sum(bkgsys[2*i+1].event_arr[bkgsys[2*i+1].event_idx]['weight']*time)
            corrbkgsysup.append(max(n1,n2)-nbkg[0])
            corrbkgsysdn.append(min(n1,n2)-nbkg[0])
            
    #compute Signal + Background wt Sys
    for s in wtsys:
        #signal
        n1 = np.sum(sigproc[0].event_arr[sigproc[0].event_idx][s+'UP']*time)
        n2 = np.sum(sigproc[0].event_arr[sigproc[0].event_idx][s+'DN']*time)
        corrsigsysup.append(max(n1,n2)-nsig[0])
        corrsigsysdn.append(min(n1,n2)-nsig[0])

        #background
        n1 = np.sum(bkgproc[0].event_arr[bkgproc[0].event_idx][s+'UP']*time)
        n2 = np.sum(bkgproc[0].event_arr[bkgproc[0].event_idx][s+'DN']*time)
        corrbkgsysup.append(max(n1,n2)-nbkg[0])
        corrbkgsysdn.append(min(n1,n2)-nbkg[0])

    #compute total events
    TOT = TOTSIG + TOTBG

    #compute mc errors
    SIGMCERR = np.sqrt(np.sum(np.array(nsig_mcerr)**2))
    BGMCERR  = np.sqrt(np.sum(np.array(nbkg_mcerr)**2))
    TOTMCERR = np.sqrt(SIGMCERR**2 + BGMCERR**2)
    
    #compute sys. errors (uncorrelated)
    SIGSYSERRUNCORR = np.sqrt(np.sum(np.array(uncorrsigsys)**2))
    BGSYSERRUNCORR  = np.sqrt(np.sum(np.array(uncorrbkgsys)**2))
    TOTSYSERRUNCORR = np.sqrt(SIGSYSERRUNCORR**2 + BGSYSERRUNCORR**2)

    #compute sys. errors (correlated)
    SIGSYSERRCORRUP = np.sqrt(np.sum(np.array(corrsigsysup)**2))
    SIGSYSERRCORRDN = np.sqrt(np.sum(np.array(corrsigsysdn)**2))
    
    BGSYSERRCORRUP  = np.sqrt(np.sum(np.array(corrbkgsysup)**2))
    BGSYSERRCORRDN  = np.sqrt(np.sum(np.array(corrbkgsysdn)**2))
    
    TOTSYSERRCORRUP = np.sqrt(SIGSYSERRCORRUP**2 + BGSYSERRCORRUP**2)
    TOTSYSERRCORRDN = np.sqrt(SIGSYSERRCORRDN**2 + BGSYSERRCORRDN**2)

    #convert sys errors to relative errors
    SIGSYSERRUNCORR /= TOTSIG
    BGSYSERRUNCORR  /= TOTBG
    TOTSYSERRUNCORR /= TOT

    SIGSYSERRCORRUP = 1 + SIGSYSERRCORRUP/TOTSIG
    SIGSYSERRCORRDN = 1 - SIGSYSERRCORRDN/TOTSIG
    
    BGSYSERRCORRUP  = 1 + BGSYSERRCORRUP/TOTBG
    BGSYSERRCORRDN  = 1 - BGSYSERRCORRDN/TOTBG
    
    TOTSYSERRCORRUP = 1 + TOTSYSERRCORRUP/TOT
    TOTSYSERRCORRDN = 1 - TOTSYSERRCORRDN/TOT

    totsum = (TOT, TOTMCERR, TOTSYSERRUNCORR, TOTSYSERRCORRUP, TOTSYSERRCORRDN)
    bkgsum = (TOTBG, BGMCERR, BGSYSERRUNCORR, BGSYSERRCORRUP, BGSYSERRCORRDN)
    sigsum = (TOTSIG, SIGMCERR, SIGSYSERRUNCORR, SIGSYSERRCORRUP, SIGSYSERRCORRDN)
    return totsum, bkgsum, sigsum

def get_DataCounts(data, cutreg):
    '''Args:
    data (analysis object): Analysis object for data
    cutreg (tuple): Cut region where we will compute the observed events
    '''
    data.reset_cut()
    data.apply_energyCut()
    data.apply_zenithCut()
    if cutreg[0]=='boxcut':
        data.apply_boxCut(cutreg[1][0], cutreg[1][1],
                        amax=cutreg[1][2],
                        bmax=cutreg[1][3])
    elif cutreg[0] in ['ellipse','hyperbola']:
        data.apply_curveCut(cutreg[1][0], cutreg[1][1],
                               cutreg[0])
    nobs = len(data.event_arr[data.event_idx])
    return nobs

#This function generates the output pickle file that can be used directly for the 
# final fit (Histfitter)
def generate_mcPickle(reglist, snom, bnom, ssys, bsys, wsys,
                     time, outfile, asimovfile=None):
    '''Args:
    reglist (dict): dictionary containing all the cur regions
    snom (list): List of analysis objects for all nominal signals
    bnom (list): List of analysis objects for all backgrounds
    ssys (list): List of all signal detector systematcis analysis objects
    bsys (list): List of all background detector systematics analsys objects
    wsys (list): List of strings for weight based systematics
    time (float): Live time for which we compute the expected out values
    outfile (str): Name of the output MC pickle file
    asimovfile (str): If a pickle filename given, the function will generate
                    all possible asimov data points and store in the given
                    output pickle file
    '''
    outdict = {key:{} for key in reglist.keys()}
    if asimovfile is not None:
        adict = {key:{} for key in reglist.keys()}
    inkeys = ['NEvents', 'MCUnc', 'UnCorr', 'CorrUp', 'CorrDn']
    for key,reg in reglist.items():
        tot, bkg, sig = get_MCCounts(snom, bnom, ssys, bsys, wsys, reg, time)
        outdict[key]['Background'] = {}
        outdict[key]['Signal'] = {}
        for i,ikey in enumerate(inkeys):
            outdict[key]['Background'][ikey] = bkg[i]
            outdict[key]['Signal'][ikey] = sig[i]
        if asimovfile is not None:
            adict[key]['asimov_s+b'] = tot[0]
            adict[key]['asimov_b'] = bkg[0]
            adict[key]['asimov_s'] = sig[0]
    #print (outdict)
    #print (adict)
    with open(outfile, 'wb') as mcpickle:
        pickle.dump(outdict, mcpickle)
    if asimovfile is not None:
        with open(asimovfile, 'wb') as aspickle:
            pickle.dump(adict, aspickle)
        return outdict, adict
    else:
        return outdict, None

#This function generates the output pickle file for data
def generate_dataPickle(reglist, data, outfile):
    '''Args:
    reglist (dict): Dictionary defining all the cut regions
    data (analysis obj): Data analysis object
    outfile (str): Name of the output pickle file
    '''
    ddict = {}
    for key,reg in reglist.items():
        val = get_DataCounts(data,reg)
        ddict[key] = val
    #print (ddict)
    with open(outfile, 'wb') as datapickle:
        pickle.dump(ddict, datapickle)
    return ddict
