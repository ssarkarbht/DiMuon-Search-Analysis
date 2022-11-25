#!/bin/python

'''
Author: Sourav Sarkar
Date: November 5, 2022
Email: ssarkar1@ualberta.ca
Description: This script processes final data and produces
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
(options, args) = parser.parse_args()

with open(options.CONFIG) as f:
        config = json.load(f)

print ("Processing Data...")
#Process the data
infiles = [config["Input"]["data_loc"]+f for f in config["Input"]["filenames"]]
ascores = [config["Input"]["score_loc"]+f for f in config["Input"]["Ascores"]]
bscores = [config["Input"]["score_loc"]+f for f in config["Input"]["Bscores"]]

data = AnalysisCutV2(infiles, ascores, bscores, None)

print ("Applying cut regions and getting event counts...")
#Get events in each cut regions
reglist = config["CutRegions"]
for key, value in reglist.items():
        print (f"---------------{key}-----------------")
        print (get_DataCounts(data, value))

#set up output directory
outdir = config["Output"]["outdir"]
if outdir is None:
        print (f"Output directory not specified in config. Will be storing output in current directory {cwd}")
        outdir = cwd+'/'

outfile = outdir + config["Output"]["outfile"]
print (f"Storing data event counts in pickle file output : {outfile}")
generate_dataPickle(reglist, data, outfile)

end_time = time.time()
duration = (end_time-start_time)/60.
print (f"Took {duration} minutes to run! GoodBye!")


