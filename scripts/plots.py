import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import mappy as mp
import matplotlib.font_manager as font_manager
import matplotlib
import os.path

def get_diff_list(list1, list2):
    ret = []
    for i in range(len(list1)):
        ret.append(list2[i]-list1[i])
    return ret


parser = argparse.ArgumentParser( description='Generate plots of before and after realignment data')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
parser.add_argument('--samples', required=True)
parser.add_argument('--read-lens', required=True)
args = parser.parse_args()

samples = args.samples.split(',')
read_lens = args.read_lens.split(',')

data = defaultdict(list)
count = 0
for rlen in read_lens:
    data['X'].append(count)
    data['XString'].append(str(rlen))
    count += 2

data['X'] = np.array(data['X'])

# Read and Plot Sniffles/CuteSV Before after Precision & Recall
pipeline_precision = defaultdict(list)
somrit_precision = defaultdict(list)
pipeline_recall = defaultdict(list)
somrit_recall = defaultdict(list)
for sample in samples:
    for rlen in read_lens:
        if not os.path.isfile(args.input_dir+"/simulation_results_"+sample+"_"+rlen+".txt"):
            continue
        with open(args.input_dir+"/simulation_results_"+sample+"_"+rlen+".txt",'r') as in_fp_tra:
            for line in in_fp_tra:
                row = line.strip().split('\t')
                pipeline_precision[rlen].append(float(row[9]))
                pipeline_recall[rlen].append(float(row[8]))
                somrit_precision[rlen].append(float(row[4]))
                somrit_recall[rlen].append(float(row[3]))



for cov in read_lens:
    data["pipeline_precision_mean"].append(np.mean(pipeline_precision[cov]))
    data["somrit_precision_mean"].append(np.mean(somrit_precision[cov]))
    data["pipeline_precision"].append(pipeline_precision[cov])
    data["somrit_precision"].append(somrit_precision[cov])
    data["pipeline_recall_mean"].append(np.mean(pipeline_recall[cov]))
    data["somrit_recall_mean"].append(np.mean(somrit_recall[cov]))
    data["pipeline_recall"].append(pipeline_recall[cov])
    data["somrit_recall"].append(somrit_recall[cov])


wd = 0.35
plt.bar(data['X'],data["pipeline_precision_mean"], color='r', width=wd, edgecolor='k', label="RTD-Control Precision", alpha=0.4)
plt.bar(data['X']+wd,data["pipeline_recall_mean"], color='g', width=wd, edgecolor='k', label="RTD-Control Sensitivity", alpha=0.4)
plt.bar(data['X']+2*wd,data["somrit_precision_mean"], color='b', width=wd, edgecolor='k', label="RTD-Diploid Precision", alpha=0.4)
plt.bar(data['X']+3*wd,data["somrit_recall_mean"], color='m', width=wd, edgecolor='k', label="RTD-Diploid Sensitivity", alpha=0.4)
plt.boxplot(data["pipeline_precision"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_precision"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["pipeline_recall"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["somrit_recall"], positions=data['X']+3*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('RTD-Control and RTD-Diploid Precision and Sensitivity', fontsize=20)
plt.xlabel('Read Length (bp)', fontsize=17)
plt.ylabel('Fraction', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_PR.png")
plt.clf()


