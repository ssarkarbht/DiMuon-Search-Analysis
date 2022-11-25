# DiMuon Search Analysis

***

> Author : Sourav Sarkar

> Email : <ssarkar1@ualberta.ca>

> Analysis Wiki : <https://wiki.icecube.wisc.edu/index.php/Search_for_Neutrino_DiMuon_Events>

***

## Analysis Overview

***

This repository contains the analysis code to perform dimuon search in IceCube up-going track sample data within reconstructed energy range \[1TeV - 50TeV]. The dimuon signal in this analysis primarily consists of $\nu_{\mu} (\bar{\nu}_{\mu})$ CCDIS interaction with charm production which eventually decays into a second high energetic muon. The background for this analysis is either (primarily) single muon tracks or misreconstructed track from other processes. In summary, the code implements the following steps:

*   Processes Monte Carlo (MC) and experimental datasets and applies event selection cuts in 4 final observable (reconstructed energy, reconstructed zenith, two dimuon classification scores) to compute expected/observed number of events in control region (CR), validation regions (VRs) and signal regions (SRs).

*   The processing also computes all systematics and limited MC statistical uncertainties for each region and samples (signal, background).

*   Validates the model by performing **background only fit** in CR and extrapolating the fit results in VRs.

*   Performs **discovery fit** by constraining both CR and SR to compute p-value with which we can exclude background-only (null) hypothesis. The alternative hypothesis for this fit is Standard Model (SM) expectation of signal+background.

*   Performs **exclusion fit** and **upper limit scan** to compute the signal strength parameter $\mu_s$ at 95% confidence level (CL) with which we can exclude the signal+background (null) hypothesis.

Details of the full analysis can be found in the analysis wiki.

## Running the Analysis Chain

***

The analysis chain can be divided into two parts. In the first part, we process the full MC or experimental dataset to get the expected or observed events in each region (and their respective uncertainties). The output from this part is saved into pickle files which are then used as input for the next step. In the second part, we use the pickle files to perform various final fits discussed above.

### Dependencies

Specific software and  data files required to run this analysis are not provided with this repository. The software and data files are located within the IceCube computing infrastructure (cobalt). Access to cobalt machines is thus required to run the full analysis chain and reproduce the results.

### Performing a Test Run

In order to run analysis tests, first step is to get the repository by running:

```
git clone <link> DiMuonAnalysis
```

Then go to the repository directory and set the required environment and library paths by running:

```
cd DiMuonAnalysis
source setup.sh
```

You should be able to run all the analysis code now. However, for convenience, all types of analysis code runs are put together into a single bash wrapper script `test/testrun.sh` where you can navigate through different part of the analysis runs. Launch the test run script:

```
cd test
./testrun.sh
```

Each individual steps/runs are discussed below:

### \[Step 1] Processing MC/Data Datasets

The scripts for running data processing are located at `process_data` and a configuration file is needed for each (MC/data) run which is located at `config`.

To run processing for MC datasets, follow:
```    
cd $ANALYSIS_DIR/process_data
python process_mc.py -c $ANALYSIS_DIR/config/mc_config.json -a
```

The option `-c` takes the config file which includes location of required data files and definition of each cut regions and additional information. The option `-a` enables the script to produce the Asimov dataset for computing expected sensitivity before unblinding.

To process burn sample data, execute:
```
cd $ANALYSIS_DIR/process_data
python process_data.py -c $ANALYSIS_DIR/config/burndata_config.json
```

If no output directory is specified in the config files, the output pickle files will be generated in the current run directory `process_data`. After running the above scripts, you should have 4 output files : `MCFullInput.pkl`, `MCAsimovInput.pkl`, `MCBurnInput.pkl`, `DataBurn.pkl`.

Note that these 4 output files are already provided with the repository in case you want to skip the processing runs (which runs for about \~1.5 hours) and only test the fit runs. Running the above commands will overwrite the outputs with same filenames and contents.

### \[Step 2] Performing Analysis Fits

Analysis fits are performed using a root based statistical analysis software [HistFitter](https://arxiv.org/pdf/1410.1280.pdf) (github repo: <https://github.com/histfitter/histfitter>). The software also provides pre-built docker container image. We will be using the container image for running our analysis fits.

(optional) Although the analysis scripts can be run directly from the dockerhub image, one can optionally pull the image and run the analysis scripts with it. If you wish to pull the software image in your local directory, run the following:
```
cd $ANALYSIS_DIR/histfitter
singularity pull docker://histfitter/histfitter:v1.0.1
```

You should now have a singularity image file (`sif` format) of the software.

Note: All the following fits are performed with a default random seed value `seed=100` to make the results reproducible. The seed value can be easily changed in the respective bash scripts.

In our analysis, we have two overlapping signal regions (SR1 and SR2). So we will be performing separate discovery and exclusion fit runs for each signal region.

1.  **Background Only Fit**

   (Option 1) If you have a pulled image of the software, run:
```
cd $ANALYSIS_DIR/histfitter
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_bkgFit.sh
```

(Option 2) If you want to directly execute the analysis script, run:
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_bkgFit.sh
```

At the end of the run, you should see an output like the following in your terminal:
```
    <INFO> YieldsTable: Will dump output dictionary:
    <INFO> YieldsTable: Fitted_err_Bkg: [39.064430238282235, 21.163238606028496, 6.907033718446899, 13.704971739356388, 11.291334274693284]
    <INFO> YieldsTable: Fitted_events_Bkg: [1536.7560547684059, 389.2182318287321, 59.77200761137556, 199.7261586501464, 91.51551107544655]
    <INFO> YieldsTable: MC_exp_err_Bkg: [500.24216528248746, 133.44776917817217, 20.614897087878376, 67.87990548185181, 32.11741307307352]
    <INFO> YieldsTable: MC_exp_events_Bkg: [1583.0449100935887, 400.9389164915122, 61.57242792064174, 205.741244053196, 94.27192154763941]
    <INFO> YieldsTable: TOTAL_FITTED_bkg_events: [1536.7560547684059, 389.2182318287321, 59.77200761137556, 199.7261586501464, 91.51551107544655]
    <INFO> YieldsTable: TOTAL_FITTED_bkg_events_err: [39.064430238282235, 21.163238606028496, 6.907033718446899, 13.704971739356388, 11.291334274693284]
    <INFO> YieldsTable: TOTAL_MC_EXP_BKG_err: [500.24216528248746, 133.44776917817217, 20.614897087878376, 67.87990548185181, 32.11741307307352]
    <INFO> YieldsTable: TOTAL_MC_EXP_BKG_events: [1583.0449100935887, 400.9389164915122, 61.57242792064174, 205.741244053196, 94.27192154763941]
    <INFO> YieldsTable: names: ['CR', 'VR1', 'VR2', 'VR3', 'VR4']
    <INFO> YieldsTable: nobs: [1537.0, 389.0, 46.0, 197.0, 93.0]
    <INFO> YieldsTable: Result written in: results/burnYield.tex
```

The run will also generate few output files in `$ANALYSIS_DIR/histfitter/results/` where the fit results are summarised in latex tables - `burnSyst.tex`, `burnYield.tex`.

2.  **Discovery Fit**
(Option 1) If you have a pulled image of the software, run:
- For SR1,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_discFit_SR1.sh
```

- For SR2,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_discFit_SR2.sh
```

(Option 2) For running the scripts directly,
- For SR1,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_discFit_SR1.sh
```

- For SR2,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_discFit_SR2.sh
```

At the end of SR1 fit, you should see a terminal output:
```
Results HypoTestCalculator_result: 
 - Null p-value = 0.0224 +/- 0.00209276
 - Significance = 2.00653 +/- 0.0392715 sigma
 - Number of Alt toys: 2500
 - Number of Null toys: 5000
 - Test statistic evaluated on data: 1.56054
 - CL_b: 0.0224 +/- 0.00209276
 - CL_s+b: 0.4564 +/- 0.00996191
 - CL_s: 20.375 +/- 1.95483
<INFO> ConfigMgrCPP: Now storing HypoTestInverterResult <hypo_discovery_Sig>
<INFO> ConfigMgrCPP: Now storing HypoTestResult <discovery_htr_Sig>
<INFO> ConfigMgrCPP: Done. Stored HypoTest(Inverter)Result and fit result in file <results/SR1_DiscFits_Output_hypotest.root>
```

At the end of SR2 fit, the output should be:
```
Results HypoTestCalculator_result: 
 - Null p-value = 0.0368 +/- 0.0018827
 - Significance = 1.78909 +/- 0.0233845 sigma
 - Number of Alt toys: 5000
 - Number of Null toys: 10000
 - Test statistic evaluated on data: 1.55911
 - CL_b: 0.0368 +/- 0.0018827
 - CL_s+b: 0.474 +/- 0.0070615
 - CL_s: 12.8804 +/- 0.686339
<INFO> ConfigMgrCPP: Now storing HypoTestInverterResult <hypo_discovery_Sig>
<INFO> ConfigMgrCPP: Now storing HypoTestResult <discovery_htr_Sig>
<INFO> ConfigMgrCPP: Done. Stored HypoTest(Inverter)Result and fit result in file <results/SR2_DiscFits_Output_hypotest.root>
```

Similar summary output files are generated in `$ANALYSIS_DIR/histfitter/results/`. For SR1 example, relevant summary files are `disc_SR1Syst.tex`, `disc_SR1Yield.tex`.

3.  **Exclusion Fit**
(Option 1) If you have a pulled image of the software, run:
- For SR1,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_exclFit_SR1.sh
```

- For SR2,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_exclFit_SR2.sh
```

(Option 2) For running the scripts directly,
- For SR1,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_exclFit_SR1.sh
```

- For SR2,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_exclFit_SR2.sh
```

At the end of SR1 fit, you should see a terminal output:
```
	CLs      = 0.310289 +/- 0.00877337
	CLb      = 0.7464 +/- 0.00870143
	CLsplusb = 0.2316 +/- 0.00596593
```

At the end of SR2 fit, you should see a terminal output:
```
	CLs      = 0.149293 +/- 0.00569904
	CLb      = 0.495 +/- 0.00707071
	CLsplusb = 0.0739 +/- 0.00261608
```

Similar summary output files are generated in `$ANALYSIS_DIR/histfitter/results/`. For SR1 example, relevant summary files are `excl_SR1Syst.tex`, `excl_SR1Yield.tex`.

4.  **Upper Limit Scan**

*[CAUTION :]* Upper limit scan can take several (4-5) hours to run and consumes a large amount of memory (~50GB). It should be run only when the time and computing resources are available.

(Option 1) If you have a pulled image of the software, run:
- For SR1,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_ulScan_SR1.sh
```

- For SR2,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs histfitter_v1.0.1.sif /DiMuonAnalysis/histfitter/scripts/run_ulScan_SR2.sh
```

(Option 2) For running the scripts directly,
- For SR1,
```

singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_ulScan_SR1.sh
```

- For SR2,
```
singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_ulScan_SR2.sh
```

### Reproducibility Check
The results of each fit can be visually checked by comparing the end part of the terminal outputs with the outputs provided above for each fit results (assuming running with same default random seed 100). In addition, md5sum hashes are provided for relevant output files for each test run in `$ANALYSIS_DIR/checksum.md5` with the list of files in `$ANALYSIS_DIR/validation_filelist.txt`.

Validity of the output file integrity can be checked by running:
```
cd $ANALYSIS_DIR
md5sum -c checksum.md5
```

**Note:** Validation of files for each tests can only be passed successfully (`Status : OK`) if the test is run to generate those files mentioned in the file list. In case, some tests are skipped, output file might not be present in the expected directory resulting in `FAILED` validation test.
