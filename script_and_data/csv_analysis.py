#Adam Naguib, 2015

import matplotlib 
import csv
import os
import matplotlib.pyplot as plt
import numpy as np 
import statistics as stats
from tabulate import tabulate
from pandas import pandas
from autoprotocol.container import Container
from autoprotocol.container_type import _CONTAINER_TYPES

mean_SD = []

CSVfiles = ([file for root, dirs, files in os.walk('/Users/adamnaguib/Desktop/data_analysis_python/script_and_data')
    for file in files
    if file.endswith('.csv') and "mean_SD" not in file])

def get_filename(file):
    for y in CSVfiles: 
        filename = (os.path.basename(file))[21:31]
    return filename

def get_meanCV(file):
    CSV_file = pandas.read_csv(file)
    expt_samples = len(CSV_file)
    DNAs = 7
    replicates = expt_samples/DNAs

    if "A13" in CSV_file["Well"].values:
        plate_map = Container(None, _CONTAINER_TYPES['384-pcr'])
    else:
        plate_map = Container(None, _CONTAINER_TYPES['96-pcr'])

    start = 0
    replicate_locs = []
    for i in range (0,DNAs-1):
        loc = [plate_map.humanize(s) for s in range(start, start + replicates)]
        replicate_locs.append(loc)
        start += replicates

    DNA_Ct = []
    for h in replicate_locs:    
        for x in h:
            Replicate_Ct_DNA = []
            data_source = open(file)    
            replicate_locations = h
            for line in data_source:
                split_line=line.split(',')
                wellID=split_line[0]
                Ct=split_line[3]
                for w in replicate_locations:
                    if w == wellID:
                        try:
                            Replicate_Ct_DNA.append(float(Ct))
                        except:
                            Replicate_Ct_DNA.append(0.0)
        DNA_Ct.append(Replicate_Ct_DNA)

    percentageCV = []
    for n in DNA_Ct:
        try:
            percentageCV.append(((stats.pstdev(n)/stats.mean(n))*100))
        except ZeroDivisionError as err:
            percentageCV.append(0.0)
    meanCV = stats.mean(percentageCV)
    for n in DNA_Ct:
        line = []
        line.append(stats.mean(n))
        line.append(stats.pstdev(n))
        mean_SD.append(line)
    writer = csv.writer(open('./output/mean_SD.csv', 'w'))
    writer.writerows(mean_SD)
    return meanCV

def graph(point_list, num_csvs, experiment_list):
    Expt  = range(0, num_csvs)
    PercentageCV = point_list
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1, axisbg='white')
    plt.subplots_adjust(bottom=0.2)
    plt.margins(0.05)
    plt.plot(Expt,PercentageCV, '--ko')   
    plt.ylim((0,5))
    plt.xticks (Expt, experiment_list, rotation=45)
    plt.title('Omni qPCR experimental consistency WC1', weight='bold',
          size=20, verticalalignment='baseline', y=1.05)
    plt.xlabel('Experiment', size=16, style='italic', labelpad=5, weight='bold')
    plt.ylabel('Experimental\n%CV', verticalalignment='top', labelpad=40,
         size=16, style='italic', x=1.05, weight='bold')
    plt.show()

points = []
experiment = []
num_csvs = 0

for f in sorted(CSVfiles):
    points.append(get_meanCV(f))
    experiment.append(get_filename(f))
    num_csvs += 1

graph(points, num_csvs, experiment)
