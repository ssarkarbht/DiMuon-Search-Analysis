#!/bin/bash

#This is a wrapper script for running various analsyis tests and 
#reproduce teh results.

#================== Data Processing ===================
echo "[STEP 1] : Data Processing ------------------"
echo "First step of the analysis chain is to process MC+Data files to extract expected event numbers for final fitting."
echo "Output of the processing is already provided in the repo. However, running this step will overwrite the output files."
echo "Do you want to proceed with data processing? (expected runtime: ~1.5HRs) Press 'y' to proceed, 'n' to skip."
echo "Press 0 to exit."
read runproc

if [[ $runproc == 'y' ]]
then
	echo "Running data processing."
	cd $ANALYSIS_DIR/process_data
	#running MC processing
	python process_mc.py -c $ANALYSIS_DIR/config/mc_config.json -a
	wait
	#running Burn Sample processing
	python process_data.py -c $ANALYSIS_DIR/config/burndata_config.json
	wait
elif [[ $runproc == 'n' ]]
then
	echo ">>>>>>>>>>>>>>>>>>>>>>>>> Skipping data processing. >>>>>>>>>>>>>>>>>>>>>>>>"
elif [[ $runproc == 0 ]]
then
	echo "Exiting test run!"
else
	echo "Invalid option. Exiting data processing."
fi

#=================== Analysis Fits ======================
#Function to display fit options
showFit () {
echo "[STEP 2] : Fitting (Counting Experiments) ------------------"
echo "We can now perform the follwoing fitting tests: "
echo "*** Burn Sample Background-Only Fit - To run it press 1, then Enter"
echo "*** Discovery Fit using Asimov Dataset (Expected Sensitivity) - To run it press 2, then Enter"
echo "*** Exclusion Fit using Asimov Dataset (Expected Sensitivity) - To run it press 3, then Enter"
echo "*** Upper Limit Scan using Asimov Dataset (Expected Sensitivity) - To run it press 4, then Enter"
echo "Press 0, then Enter to exit."
read testtype
if [[ $testtype -gt 1 ]] && [[ $testtype -lt 5 ]]
then
        echo "*** Press 1 for fitting Signal Region 1 (SR1); then Enter"
        echo "*** Press 2 for fitting Signal Region 2 (SR2); then Enter"
        read srtype
        if [[ $srtype =~ 1 ]] && [[ $srtype =~ 2 ]]
        then
                echo "Invalid option. Exiting test run."
                return
        fi
fi
}

#Function to run different fits
runFit () {
if [[ $1 == 1 ]]
then
	echo "Running Background-Only Burn Sample Fit"
	singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_bkgFit.sh
	wait
	echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
	read doMore
elif [[ $1 == 2 ]]
then
	if [[ $2 == 1 ]]
	then
		echo "Running Discovery (Asimov) Fit for Signal Region 1 (SR1)"
		singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_discFit_SR1.sh
		wait
		echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
		read doMore
	elif [[ $2 == 2 ]]
	then
		echo "Running Discovery (Asimov) Fit for Signal Region 2 (SR2)"
		singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_discFit_SR2.sh
		wait
                echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
                read doMore
	else
		echo "Invalid option. Exiting test run."
		return
	fi
elif [[ $1 == 3 ]]
then
        if [[ $2 == 1 ]] 
        then
		echo "Running Exclusion (Asimov) Fit for Signal Region 1 (SR1)"
		singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_exclFit_SR1.sh
		wait
                echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
                read doMore
        elif [[ $2 == 2 ]]
	then
                echo "Running Exclusion (Asimov) Fit for Signal Region 2 (SR2)"
		singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_exclFit_SR2.sh
		wait
                echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
                read doMore
	else
		echo "Invalid option. Exiting test run."
		return
        fi
elif [[ $1 == 4 ]]
then
	echo "Upper limit scan can take several (4-5) hours and consumes large memory."
	echo "Do you still want to proceed? Press 'y' to proceed, 'n' to exit."
	read dcheck
	if [[ $2 == 1 ]] && [[ $dcheck == 'y' ]]
	then
		echo "Running Upper Limit scan for Signal Region 1 (SR1)"
		singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_ulScan_SR1.sh
		wait
                echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
                read doMore
        elif [[ $2 == 2 ]] && [[ $dcheck == 'y' ]]
        then
                echo "Running Upper Limit scan for Signal Region 2 (SR2)"
		singularity exec -B $ANALYSIS_DIR:/DiMuonAnalysis -B /cvmfs:/cvmfs docker://histfitter/histfitter:v1.0.1 /DiMuonAnalysis/histfitter/scripts/run_ulScan_SR2.sh
		wait
                echo "Completed current fit. To continue with more tests, press 'c'. To exit, press 0"
                read doMore
	else
		echo "Invalid option. Exiting test run."
		return
	fi
elif [[ $1 == 0 ]]
then
	echo "Exiting test run!"
	return
else
	echo "Invalid option. Exiting Fitting."
	return
fi
}

if [[ $runproc != 0 ]]
then
	showFit
	runFit $testtype $srtype
fi

while [[ $doMore == 'c' ]] && [[ $testtype != 0 ]]
do
	showFit
	runFit $testtype $srtype
done

echo "Current Testing Session Compteled."
