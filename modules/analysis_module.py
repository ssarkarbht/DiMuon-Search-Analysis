#!/bin/python

'''
Author: Sourav Sarkar
Email: ssarkar1@ulaberta.ca
Date: 20 August, 2022
Description: This script contains all the final analysis functions,
	and modules to update event weights, make analysis cuts and
	distributions, and get final level results.
'''

#import stuff
import nuSQUIDSpy as nsq
import numpy as np
import h5py as h5
from scipy.interpolate import interp1d

#convert pdg code to nuSQuIDS compatible particle id
def convert_pdg_nsq(pdg):
	if pdg==12:#nuE
		return (0,0)
	elif pdg==-12:#nuEBar
		return (0,1)
	elif pdg==14:#nuMu
		return (1,0)
	elif pdg==-14:#nuMuBar
		return (1,1)
	elif pdg==16:#nuTau
		return (2,0)
	elif pdg==-16:#nuTauBar
		return (2,1)
	else:
		print (f"Unkown neutrino pdg type {pdg}")
		assert False
	return

#convert pdg code to nuflux compatible particle id
def convert_pdg_nuflux(pdg):
	if pdg==12:#nuE
		return nuflux.NuE
	elif pdg==-12:#nuEBar
		return nuflux.NuEBar
	elif pdg==14:#nuMu
		return nuflux.NuMu
	elif pdg==-14:#nuMuBar
		return nuflux.NuMuBar
	elif pdg==16:#nuTau
		return nuflux.NuTau
	elif pdg==-16:#nuTauBar
		return nuflux.NuTauBar
	else:
		print (f"Unkown neutrino pdg type {pdg}")
		assert False
	return

#Cross section module
class CrossSection:
	''' This class handles all the final cross-sections
	from both single muon and dimuon processes.
	'''
	def __init__(self, xsdir, ct18=True):
		'''
		Args: xsdir (str): Parent directory that contains
			all cross-section data
		      ct18 (bool): If we want to make corrections
			in the final cross-section from CT14nlo
			charm fraction to CT18Anlo cross-section.
			default=True.
		'''
		#get data dir of different processes
		#CCDIS
		self.ccdis_path = xsdir+'CSMS/'
		#Neutrino Trident Cross-sections
		self.ntp_path   = xsdir+'ZhouBeacom/'
		#Charm fraction and Charm Muon cross-sections
		self.charm_path = xsdir+'CharmFraction/'

		#build the ccdis interpolation
		en_points = np.linspace(1,12,111)
		nu_cc_iso = np.loadtxt(self.ccdis_path+\
				'total_nu_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',
				skiprows=1)
		nubar_cc_iso = np.loadtxt(self.ccdis_path+\
				'total_nubar_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',
				skiprows=1)

		self.nu_cc_func = interp1d(en_points, nu_cc_iso, kind='cubic')
		self.nubar_cc_func = interp1d(en_points, nubar_cc_iso, kind='cubic')

		#build the trident interpolation
		elastic   = np.loadtxt(self.ntp_path+'NuMu_CCNC_Elastic_O16.dat',
						delimiter=',', comments='#')
		disquark  = np.loadtxt(self.ntp_path+'NuMu_CCNC_DISQuark_O16.dat',
						delimiter=',', comments='#')
		disphoton = np.loadtxt(self.ntp_path+'NuMu_CCNC_DISPhoton_O16.dat',
						delimiter=',', comments='#')

		elastic_logen = np.log10(elastic[:,0])
		elastic_logxs = np.log10(elastic[:,1])
		self.elastic_func = interp1d(elastic_logen, elastic_logxs,
								kind='slinear')

		disquark_logen = np.log10(disquark[:,0])
		disquark_logxs = np.log10(disquark[:,1])
		self.disquark_func = interp1d(disquark_logen, disquark_logxs,
								kind='slinear')

		disphoton_logen = np.log10(disphoton[:,0])
		disphoton_logxs = np.log10(disphoton[:,1])
		self.disphoton_func = interp1d(disphoton_logen, disphoton_logxs,
								kind='slinear')

		#build the charm muon fraction interpolation
		numu_charmfrac = np.loadtxt(self.charm_path+'NuMu_CharmMuon_CC_Fraction.dat',
								comments='#')
		logen     = np.log10(numu_charmfrac[:,0])
		nufrac    = numu_charmfrac[:,1]
		nubarfrac = numu_charmfrac[:,2]
		if ct18:#make correction to charm fraction (ct14 to ct18)
			ct18data = np.loadtxt(self.charm_path+'CT18A_Nu+NuBar_CC_ISO_Charm_Fraction.dat',
							delimiter=',', comments='#')
			tempen = np.log10(ct18data[:,0])
			nu_ratio = ct18data[:,3]
			nubar_ratio = ct18data[:,4]
			self.nu_frac_ratio = interp1d(tempen, nu_ratio, kind='cubic')
			self.nubar_frac_ratio = interp1d(tempen, nubar_ratio, kind='cubic')

			nufrac *= self.nu_frac_ratio(logen)
			nubarfrac *= self.nubar_frac_ratio(logen)

		self.nu_charmmuon_frac = interp1d(logen, nufrac, kind='slinear',
								bounds_error=False, fill_value='extrapolate')
		self.nubar_charmmuon_frac = interp1d(logen, nubarfrac, kind='slinear',
								bounds_error=False, fill_value='extrapolate')
		
	def totalCCDIS(self, nuen, nupdg):
		'''Computes the total CCDIS xs for a given energy and pdg
		on isoscalar target nucleon.
		'''
		logen = np.log10(nuen)
		xsval = np.zeros(len(nuen))
		nuidx = np.where(nupdg>0)[0]
		nubaridx = np.where(nupdg<0)[0]

		#nu
		xsval[nuidx] = self.nu_cc_func(logen[nuidx])*1e9 #pb (1 pb=10^-36 cm^2)
		#nubar
		xsval[nubaridx] = self.nubar_cc_func(logen[nubaridx])*1e9
		assert len(nuidx)+len(nubaridx)==len(xsval), "Mismatch in index!"
		return xsval


	def totalNTPCCNC(self, nuen):#cross-section is same for nu/nubar
		'''Computes the total trident cross-section from numu/numubar
		cc+nc channel on isoscalar target nucleon
		'''
		logen = np.log10(nuen)
		elasticxs   = 10**self.elastic_func(logen)
		disquarkxs  = 10**self.disquark_func(logen)
		disphotonxs = 10**self.disphoton_func(logen)

		totalxs = elasticxs+disquarkxs+disphotonxs
		iso_xs = totalxs*nuen/16.*1e36 #pb
		return iso_xs

	def CharmMuonFraction(self, nuen, nupdg):
		'''Computes the fractional probability of producing muon
		from outgoing charm quark in ccdis interaction.
		'''
		logen = np.log10(nuen)
		fracval = np.zeros(len(nuen))
		nuidx = np.where(nupdg>0)[0]
		nubaridx = np.where(nupdg<0)[0]

		fracval[nuidx] = self.nu_charmmuon_frac(logen[nuidx])
		fracval[nubaridx] = self.nubar_charmmuon_frac(logen[nubaridx])
		assert len(nuidx)+len(nubaridx)==len(fracval), "Mismatch in index!"
		return fracval

	def extract_charmFraction(self, xs_in, nuen, nupdg):
		'''This function extracts charm muon fraction for each simulated 
		charm dimuon events for the given total cross-section array
		xs_in (array) : final cross-section in pb
		'''
		ccxs = self.totalCCDIS(nuen, nupdg)
		xs_in *= 1e36 #convert to pb
		return xs_in/ccxs

	def get_pdfset_correction(self, nuen, nupdg):
		''' Computes a correction factor to the cross-section for
		updating the pdfset from CT14nlo to CT18Anlo
		'''
		logen = np.log10(nuen)
		corr_factor = np.ones(len(nuen))
		nuidx = np.where(nupdg>0)[0]#nu
		nubaridx = np.where(nupdg<0)[0]#nubar

		corr_factor[nuidx] = self.nu_frac_ratio(logen[nuidx])
		corr_factor[nubaridx] = self.nubar_frac_ratio(logen[nubaridx])
		assert len(nuidx)+len(nubaridx)==len(corr_factor), "Mismatch in index!"
		return corr_factor
		
class Reweighting:
	'''
	This class takes the one weight arrays for a dataset and
	makes correction/updates the definition to match across
	all datasets. Also reweights according to the flux and
	cross section systematics.
	'''
	def __init__(self, filename, xsdir, nsqdir):
		'''
		Args: filename (str): Name of the input MC event file (HDF5 format)
		      xsdir (str): Parent directory containing cross-section data
		      nsqdir (str): Directory name containing nuSQuIDs files
		'''
		#load the hdf5 file
		hfile = h5.File(filename,"r")
		#get event ID
		self.event_id = hfile["EventID"][:]
		#get run ID
		self.run_id   = hfile["RunID"][:]
		#get Neutrino energy, zenith and PDG code
		self.nuen  = hfile["NeutrinoEnergy"][:]
		self.nuzen = hfile["NeutrinoZenith"][:]
		self.nupdg = hfile["NeutrinoPDGCode"][:]
		#get Weight information
		self.cdepth= hfile["TotalColumnDepth"][:]
		self.onew  = hfile["OneWeight"][:]
		#Build the Cross-section object
		self.cross_section  = CrossSection(xsdir)
		#nuSQuIDs directory
		self.nsq_dir = nsqdir

	def interaction_weight(self, xs, convert_unit=False):
		'''Computes the interaction weight of an event.
		Args: xs (arr): cross-section in pb or cm^2
		      convert_unit (bool): if cross-section input
		is in pb, we need to convert it in cm^2. Default=False
		'''
		#Constants
		NA = 6.022140857e+23 #Avogadro's number
		conversion = 1e-36 #pb to cm^2
		if convert_unit:
			exp_factor = -1.*self.cdepth*NA*xs*conversion
		else:
			exp_factor = -1.*self.cdepth*NA*xs
		val = 1-np.exp(exp_factor)
		#if cross-section is too small, it can result to zero
		#interaction weight due to rounding error. In this case,
		#we switch to approximate value
		zero_idx = np.where(val==0.)[0]
		if len(zero_idx)>0:
			val[zero_idx] = -1.*exp_factor[zero_idx]
		return val

	def update_interaction_weight(self, xs):
		'''This function implements the correction to the interaction
		 weight that was following
		NuGen method, however, totalcolumndepth definition is not the
		same for NuGen and LI, thus need weight correction to some of
		the charm dataset
		'''
		#Constants
		ice_density = 0.92 #g/cm^3 denisty of ice 
		proton_mass = 1.6726219e-24 #gm
		NA          = 6.022140857e+23 #Avogadro's number		
		#Compute the old (incorrect) interaction weight factor
		exp_factor_old = -1.*self.cdepth*ice_density*xs/proton_mass
		val_old = 1-np.exp(exp_factor_old)
		zero_idx = np.where(val_old==0.)[0]
		if len(zero_idx)>0:
			val_old[zero_idx] = -1.*exp_factor_old[zero_idx]
		#Now compute the new (correct) interaction weight factor
		val_new = self.interaction_weight(xs)
		#replace the old factor with new one in one weight value
		correction_factor = val_new/val_old
		self.onew *= correction_factor
		return correction_factor

	def charm_pdf_correction(self, xs):
		'''This function performs the pdfset correction from CT14nlo
		to CT18Anlo for a given corss-section array (for Charm Muon MC)
		'''
		val_old = self.interaction_weight(xs)
		newxs = xs*self.cross_section.get_pdfset_correction(self.nuen, self.nupdg)
		val_new = self.interaction_weight(newxs)

		correction_factor = val_new/val_old
		self.onew *= correction_factor
		return correction_factor

	def charm_xs_sys(self, xs, sys):
		'''This function computs the systematics variation due to
		cross-section uncertainity (+-10%)
		'''
		val_old = self.interaction_weight(xs)
		if sys=='up':
			val_new = self.interaction_weight(xs*1.1)
		if sys=='dn':
			val_new = self.interaction_weight(xs*0.9)

		correction_factor = val_new/val_old
		self.onew *= correction_factor
		return correction_factor

	def load_xsdata(self, xsfile=None, runnum=None):
		'''Get the total cross-section data for each individual
		events (Charm MC)
		'''
		# if loading xs for CCDIS compute it internally
		if (xsfile is None) and (runnum is None):
			return self.cross_section.totalCCDIS(self.nuen, self.nupdg)
		# if loading xs for charm muons, get it from external files
		hfile = h5.File(xsfile,"r")
		xsdata = hfile[runnum][:]
		xsidx = xsdata['idx'].astype(str)
		xsval = xsdata['xs']

		datadict = dict(zip(xsidx,xsval))
		xsarr = []
		for idx in range(len(self.event_id)):
			tag = str(self.run_id[idx])+\
				'{:05d}'.format(self.event_id[idx])
			xsarr.append(datadict[tag])
		return np.array(xsarr)

	def undo_propagation_weight(self):
		'''Initially the weight module implemented the propagation
		weight separatley from flux models. Since the MEOWS dataset
		has nuflux and nusquids propagation coupled together, we need
		to decouple our fluxless propagation weight.
		'''
		#old fluxless nuSQuIDs files to be used in decoupling
		numu_nsq = nsq.nuSQUIDSAtm(self.nsq_dir+\
				"numu_nsq_propagation_weight_latest.h5")
		nutau_nsq = nsq.nuSQUIDSAtm(self.nsq_dir+\
				"nutau_nsq_propagation_weight_latest.h5")
		units = nsq.Const()
		for idx in range(len(self.nuen)):
			inistate = 1.0e18*(self.nuen[idx]*units.GeV)**(-1.)
			nuflav, nutype = convert_pdg_nsq(self.nupdg[idx])
			if self.nuen[idx]>1e8: self.nuen[idx] = 1e8 #GeV
			if nuflav==1:
				propw = numu_nsq.EvalFlavor(nuflav, np.cos(self.nuzen[idx]),
						self.nuen[idx]*units.GeV, nutype)
			elif nuflav==2:
				propw = nutau_nsq.EvalFlavor(nuflav, np.cos(self.nuzen[idx]),
						self.nuen[idx]*units.GeV, nutype)
			else:
				print ("Shouldn't have other nu flavour correction!")
				assert False
			propw /= inistate
			self.onew[idx] /= propw
			
		return

	def update_cross_section(self, sys=None):
		''' CC DIS total cross section is inclusive of contribution
		from all quarks, therefore total cross section from CSMS in this
		process leads to double counting when the charm muon dataset
		is included in the analysis. We thus need to update the CCDIS
		single muon dataset cross-sections. In addition, this function
		computes the systematics variation due to cross-section uncertainties
		(+-10%) for CCDIS single muon MC sets.
		'''
		xs_original = self.cross_section.totalCCDIS(self.nuen, self.nupdg)
		intw_original = self.interaction_weight(xs_original, convert_unit=True)
		if sys is None:
			xs_update = xs_original*(1-self.cross_section.CharmMuonFraction(self.nuen,
										self.nupdg))
		if sys=='up':
			xs_update = 1.1*xs_original*(1-self.cross_section.CharmMuonFraction(self.nuen,
					self.nupdg))
		if sys=='dn':
			xs_update = 0.9*xs_original*(1-self.cross_section.CharmMuonFraction(self.nuen,
					self.nupdg))

		intw_update = self.interaction_weight(xs_update, convert_unit=True)

		update_factor = intw_update/intw_original
		self.onew *= update_factor
		return

	def get_conv_flux(self, flux='Atmospheric_Conventional.h5'):
		'''Apply atmospheric conventional flux model
		Default : H3A+Sibyll2.3c
		'''
		nsq_flux = nsq.nuSQUIDSAtm(self.nsq_dir+flux)
		units = nsq.Const()
		for i in range(len(self.nuen)):
			nuflav, nutype = convert_pdg_nsq(self.nupdg[i])
			fluxw = nsq_flux.EvalFlavor(nuflav, np.cos(self.nuzen[i]),
								self.nuen[i]*units.GeV, nutype)
			self.onew[i] *= fluxw
		return

	def get_prompt_flux(self, flux='Atmospheric_Prompt.h5'):
		'''Apply atmospheric prompt flux model
		Default: BERSS
		'''
		nsq_flux = nsq.nuSQUIDSAtm(self.nsq_dir+flux)
		units = nsq.Const()
		for i in range(len(self.nuen)):
			nuflav, nutype = convert_pdg_nsq(self.nupdg[i])
			fluxw = nsq_flux.EvalFlavor(nuflav, np.cos(self.nuzen[i]),
								self.nuen[i]*units.GeV, nutype)
			self.onew[i] *= fluxw
		return

	def get_astro_flux(self, flux='Astrophysical.h5'):
		'''Apply astrophysical flux model
		Default: Astro NuMu 9.5 years fit results
		'''
		nsq_flux = nsq.nuSQUIDSAtm(self.nsq_dir+flux)
		units = nsq.Const()
		for i in range(len(self.nuen)):
			nuflav, nutype = convert_pdg_nsq(self.nupdg[i])
			fluxw = nsq_flux.EvalFlavor(nuflav, np.cos(self.nuzen[i]),
								self.nuen[i]*units.GeV, nutype)
			self.onew[i] *= fluxw
		return

	def get_total_flux(self, convflux='Atmospheric_Conventional.h5',
			promptflux='Atmospheric_Prompt.h5',
			astroflux='Astrophysical.h5',
			conv_norm=None, conv_gamma=None,
			astro_norm=None, astro_gamma=None):
		'''Apply total flux. In addtion, this function also applies
		systematics variation (flux normalization + spectral index)
		for atmospheric convetional and astrophysical flux models.
		'''
		flux1 = nsq.nuSQUIDSAtm(self.nsq_dir+convflux)
		flux2 = nsq.nuSQUIDSAtm(self.nsq_dir+promptflux)
		flux3 = nsq.nuSQUIDSAtm(self.nsq_dir+astroflux)
		units = nsq.Const()
		for i in range(len(self.nuen)):
			nuflav, nutype = convert_pdg_nsq(self.nupdg[i])
			flux1w = flux1.EvalFlavor(nuflav, np.cos(self.nuzen[i]),
								self.nuen[i]*units.GeV, nutype)
			if conv_norm=='up':
				flux1w *= 1.3
			if conv_norm=='dn':
				flux1w *= 0.7
			if conv_gamma=='up':
				flux1w *= (self.nuen[i]/1e3)**(0.01)
			if conv_gamma=='dn':
				flux1w *= (self.nuen[i]/1e3)**(-0.01)

			flux2w = flux2.EvalFlavor(nuflav, np.cos(self.nuzen[i]),
								self.nuen[i]*units.GeV, nutype)

			flux3w = flux3.EvalFlavor(nuflav, np.cos(self.nuzen[i]),
								self.nuen[i]*units.GeV, nutype)
			if astro_norm=='up':
				flux3w *= 1.18
			if astro_norm=='dn':
				flux3w *= 0.82
			if astro_gamma=='up':
				flux3w *= (self.nuen[i]/1e5)**(0.09)
			if astro_gamma=='dn':
				flux3w *= (self.nuen[i]/1e5)**(-0.09)

			total_flux = flux1w+flux2w+flux3w
			self.onew[i] *= total_flux
		return



class AnalysisCut:
	def __init__(self, file):
		'''Parameters:
		filename (str) : name of the file to process
		'''
		#get the h5 file
		self.file = h5.File(file, "r")
		#get mask for event idxs
		eprop = np.array(self.file['RunID'])
		self.event_idx = np.ones(len(eprop), dtype=np.bool)
		#get the mask for DOM idxs
		dprop = np.array(self.file['ResidualTime'])
		self.dom_idx = np.ones(len(dprop), dtype=np.bool)

		#get event to dom idxs mapping
		vnum = np.array(self.file["VertexNumber"])
		self.end_idx = np.cumsum(vnum)
		self.start_idx = self.end_idx-vnum
		
	def apply_energyCut(self, emin=300.0, emax=3e4):
		ereco = np.array(self.file['RecoEnergy'])
		#get the event idxs to cut
		idx_removed = np.where(np.logical_or((ereco<emin), (ereco>emax)))[0]
		self.event_idx[idx_removed] = False
		#cut the DOM level idxs
		for i in idx_removed:
			self.dom_idx[self.start_idx[i]:self.end_idx[i]] = False
		return None

	def apply_zenithCut(self, zmin=-1.0, zmax=0.0):
		zreco = np.cos(np.array(self.file['RecoZenith']))
		idx_removed = np.where(np.logical_or((zreco<zmin), (zreco>zmax)))[0]
		self.event_idx[idx_removed] = False
		for i in idx_removed:
			self.dom_idx[self.start_idx[i]:self.end_idx[i]] = False
		return None

	def apply_ndomCut(self, nmin=10, nmax=500):
		ndoms = np.array(self.file['VertexNumber'])
		idx_removed = np.where(np.logical_or((ndoms<nmin), (ndoms>nmax)))[0]
		self.event_idx[idx_removed] = False
		for i in idx_removed:
			self.dom_idx[self.start_idx[i]:self.end_idx[i]] = False
		return None


class AnalysisCutV2:
	'''This class handles the final analysis cuts to get expected/observed
	events in each cut regions. The cuts implemented in this class includes
	reco. energy, reco. zenith, GNetA score, GNetB score. The class also
	handles all the weight based systematics variation in a coherent way
	(i.e. avoiding multiple copies of the same nominal array information)
	'''
	def __init__(self, files, Ascores, Bscores, weights):
		'''
		Args: files (list): list of files containg MC events
		      Ascores (list): list of files containing GNetA scores
		      BScores (list): List of files containing GNetB scores
		      weights (arr): Final weight array computed using
			reweighting class.
		'''
		print ("Processing files for building final arrays... Will take a while.")
		self.filenames = files
		#get reconstructed energy, zenith and class label
		recoTE = np.array([])
		recoCZ = np.array([])
		for f in files:
			hf = h5.File(f,'r')
			e = np.array(hf['RecoEnergy'])
			z = np.cos(np.array(hf['RecoZenith']))

			recoTE = np.append(recoTE, e)
			recoCZ = np.append(recoCZ, z)
		#if there is only one score file, load the score directly
		if len(Ascores)>1:
			#get A scores from the files
			score_A = self._merge_score(self._load_score(Ascores))
			#get B scores from the files
			score_B = self._merge_score(self._load_score(Bscores))
		#otherwise merge the scores from multple files
		else:
			score_A = self._load_score(Ascores)[0]
			score_B = self._load_score(Bscores)[0]

		#get mapping for both scores
		idxa, idx = self._combine_score(score_A, score_B)

		#total events
		self.NEvents = len(recoTE)
		#initialize event mask array
		self.event_idx = np.ones(self.NEvents, dtype=np.bool)
		#the following mask is only for EZ cuts
		self.EZevent_idx = np.ones(self.NEvents, dtype=np.bool)
		#build the final analysis array
		self.event_arr = np.zeros(self.NEvents,
			dtype=[('id',int), ('class',int), ('weight', float),
					('apred',float), ('bpred',float),
					('etrack',float), ('cosz',float)])
		#fill the array
		self.event_arr['id'] = idxa
		self.event_arr['class'] = score_A['class']
		self.event_arr['weight'] = weights
		self.event_arr['apred'] = score_A['prediction']
		self.event_arr['bpred'] = score_B['prediction']
		self.event_arr['etrack'] = recoTE
		self.event_arr['cosz'] = recoCZ
		print ("Done building arrays!")

	def _load_score(self, files):
		'''Load the scores from each file and put it in a list
		'''
		holder = []
		for file in files:
			arr = np.loadtxt(file, delimiter=',',
				dtype=[('runID', int), ('eventID', int), ('class', int), ('prediction', float)])
			holder.append(arr)
		return holder

	def _merge_score(self, scorelist):
		'''Merge the list of score arrays 
		'''
		scores = tuple(scorelist)
		return np.concatenate(scores,axis=0)

	def _combine_score(self, score1, score2):
		'''Build unique ids for each event from run and event id
		and use it to map GNetA and GNetB score for each event
		'''
		idxa = np.char.add(score1['runID'].astype('str'),
					score1['eventID'].astype('str')).astype(int)
		idxb = np.char.add(score2['runID'].astype('str'),
					score2['eventID'].astype('str')).astype(int)
		idx = np.where(np.in1d(idxb, idxa))[0]
		return idxa, idx

	def apply_energyCut(self, emin=1.1e3, emax=5e4):
		'''Apply final analysis reconstruicted energy cuts
		'''
		ereco = self.event_arr['etrack']
		#get the event idxs to cut
		idx_removed = np.where(np.logical_or((ereco<emin), (ereco>emax)))[0]
		self.event_idx[idx_removed] = False
		self.EZevent_idx[idx_removed] = False

		#cut the DOM level idxs
		#for i in idx_removed:
		#	self.dom_idx[self.start_idx[i]:self.end_idx[i]] = False
		return None

	def apply_zenithCut(self, zmin=-1.0, zmax=0.0):
		'''Apply final analysis recosntructed zenith cut
		'''
		zreco = self.event_arr['cosz']
		idx_removed = np.where(np.logical_or((zreco<zmin), (zreco>zmax)))[0]
		self.event_idx[idx_removed] = False
		self.EZevent_idx[idx_removed] = False

		#for i in idx_removed:
		#	self.dom_idx[self.start_idx[i]:self.end_idx[i]] = False
		return None

	def apply_boxCut(self, amin, bmin, amax=1.0, bmax=1.0):
		'''Apply a box cut on the 2D GNetA and GNetB score
		observable space.
		'''
		ascore = self.event_arr['apred']
		bscore = self.event_arr['bpred']
		#cut on a range
		idx_removed = np.where(np.logical_or((ascore<amin), (ascore>amax)))[0]
		self.event_idx[idx_removed] = False

		#cut on b range
		idx_removed = np.where(np.logical_or((bscore<bmin), (bscore>bmax)))[0]
		self.event_idx[idx_removed] = False

		return None

	def _ellipse(self,a,b,alpha):
		'''Elliptical function that computes the contour line
		value for give GNetA and GNetB scores (a,b) and the shape
		factor alpha for the elliptical curve.
		'''
		#range normalization factors (based on ROI range)
		SA = 5.0#1/(1-0.8)
		SB = 1.67#1/(1-0.4)
		val = 1 - SA*alpha*(a-1)**2 - SB*(1-alpha)*(b-1)**2
		return val

	def _hyperbola(self,a,b,alpha,beta=20):
		'''Hyperbolic function that computes the contour line value
		for give GNetA and GNetB scores (a,b) and the shape factors
		alpha and beta. Default value of beta = 20.0
		'''
		SA = 5.0
		SB = 1.67
		Norm = (SA-SB)*alpha+SB
		val = (SA*alpha*a+SB*(1-alpha)*b-beta*(a-1)*(b-1))/Norm
		return val

	def apply_curveCut(self, tmin, alpha, curve):
		'''Apply curve cut (either elliptical or hyperbolic) at a
		give curve value (tmin), curve shape factor (alpha).
		curve (str) : 'ellipse' or 'hyperbola'
		'''
		ascore = self.event_arr['apred']
		bscore = self.event_arr['bpred']
		#compute t values
		if curve=='ellipse':
			tval = self._ellipse(ascore, bscore, alpha)
		elif curve=='hyperbola':
			tval = self._hyperbola(ascore, bscore, alpha)
		#cut on t value
		idx_removed = np.where(tval<tmin)[0]
		self.event_idx[idx_removed] = False
		return None

	def reset_cut(self, keepEZCut=False):
		'''Reset all the existing cut, if
		keepEZCut (bool) is True only score cuts are removed.
		'''
		if keepEZCut:
			self.event_idx = self.EZevent_idx
		else:
			self.event_idx = np.ones(self.NEvents, dtype=np.bool)
		return None

	def get_subclasses(self):
		'''For dimuon events, get idxs of each classes (A,B and C)
		'''
		seglen = np.array([])
		seplen = np.array([])
		for f in self.filenames:
			hf = h5.File(f,"r")
			s1 = hf['SegLen'][::2]
			s2 = hf['SegLen'][1::2]
			seg = np.minimum(s1,s2)
			sep = hf['SepLen'][1::2]
        
			seglen = np.append(seglen, seg)
			seplen = np.append(seplen, sep)

		self.a_events = np.zeros(self.NEvents, dtype=np.bool)
		idxa = np.where((seglen>=200.)&(seplen>=25.))[0]
		self.a_events[idxa] = True

		self.b_events = np.zeros(self.NEvents, dtype=np.bool)
		idxb = np.where(seglen>=200.)[0]
		self.b_events[idxb] = True

		self.c_events = np.zeros(self.NEvents, dtype=np.bool)
		idxc = np.where(seglen<200.)[0]
		self.c_events[idxc] = True

		return None

	def get_subclass_cut(self):
		'''Get the events in each individual classes after
		applying all the analysis cuts
		'''
		cutidx = np.where(self.event_idx)
		aidx = np.where(self.a_events)
		bidx = np.where(self.b_events)
		cidx = np.where(self.c_events)

		acut = np.intersect1d(aidx, cutidx)
		bcut = np.intersect1d(bidx, cutidx)
		ccut = np.intersect1d(cidx, cutidx)
		return acut, bcut, ccut

	def add_systematics_weight(self, arr, systype):
		'''This function adds the systematics weights to the 
		initial array.
		'''
		old_dtype = self.event_arr.dtype.descr
		new_dtype = np.dtype(old_dtype + [(systype, float)])
		dummy_arr = np.zeros(self.event_arr.shape, dtype=new_dtype)
		for key in old_dtype:
			dummy_arr[key[0]] = self.event_arr[key[0]]

		dummy_arr[systype] = arr
		self.event_arr = dummy_arr
		return None
