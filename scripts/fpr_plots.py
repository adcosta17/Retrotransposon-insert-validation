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
parser.add_argument('--input', required=True)
parser.add_argument('--data', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
args = parser.parse_args()


data = defaultdict(list)
count = 0
splits = []
for i in range(3):
    data['X'].append(count)
    text = "HPRC-"
    splits_text = "HPRC-"
    if i == 0:
        text += "Long (29.1 kbp)"
        splits_text += "Long"
    elif i == 1:
        text += "Medium (9.7 kbp)"
        splits_text += "Medium"
    else:
        text += "Short (6.0 kbp)"
        splits_text += "Short"
    splits.append(splits_text)
    data['XString'].append(text)
    count += 6

print(splits)
data['X'] = np.array(data['X'])

# Read and Plot Sniffles/CuteSV Before after Precision & Recall
before_line_rate = defaultdict(list)
after_line_rate = defaultdict(list)
before_alu_rate = defaultdict(list)
after_alu_rate = defaultdict(list)
before_sva_rate = defaultdict(list)
after_sva_rate = defaultdict(list)
with open(args.input, 'r') as in_data:
    count = 0
    for line in in_data:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        #print(row)
        if row[2] == "old":
            before_line_rate[row[1]].append(float(row[4]))
            before_alu_rate[row[1]].append(float(row[6]))
            before_sva_rate[row[1]].append(float(row[8]))
        else:
            after_line_rate[row[1]].append(float(row[4]))
            after_alu_rate[row[1]].append(float(row[6]))
            after_sva_rate[row[1]].append(float(row[8]))
        

test_coverages = defaultdict(list)
with open(args.data, 'r') as in_data:
    count = 0
    for line in in_data:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        #print(row)
        test_coverages[row[1]].append(int(row[7]))

for cov in splits:
    data["before_line_mean"].append(np.mean(before_line_rate[cov]))
    data["before_line"].append(before_line_rate[cov])
    data["before_alu_mean"].append(np.mean(before_alu_rate[cov]))
    data["before_alu"].append(before_alu_rate[cov])
    data["before_sva_mean"].append(np.mean(before_sva_rate[cov]))
    data["before_sva"].append(before_sva_rate[cov])
    data["after_line_mean"].append(np.mean(after_line_rate[cov]))
    data["after_line"].append(after_line_rate[cov])
    data["after_alu_mean"].append(np.mean(after_alu_rate[cov]))
    data["after_alu"].append(after_alu_rate[cov])
    data["after_sva_mean"].append(np.mean(after_sva_rate[cov]))
    data["after_sva"].append(after_sva_rate[cov])
    data["test_coverage"].append(test_coverages[cov])



wd = 1
fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
ax1.bar(data['X'],data["before_line_mean"], color='r', width=wd, edgecolor='k', label="RTD-Control LINE FPR", alpha=0.4)
ax1.bar(data['X']+wd,data["before_alu_mean"], color='g', width=wd, edgecolor='k', label="RTD-Control Alu FPR", alpha=0.4)
ax1.bar(data['X']+2*wd+0.5,data["after_line_mean"], color='b', width=wd, edgecolor='k', label="RTD-Diploid LINE FPR", alpha=0.4)
ax1.bar(data['X']+3*wd+0.5,data["after_alu_mean"], color='m', width=wd, edgecolor='k', label="RTD-Diploid Alu FPR", alpha=0.4)
ax1.boxplot(data["before_line"], positions=data['X'], widths=wd)
ax1.boxplot(data["before_alu"], positions=data['X']+wd, widths=wd)
ax1.boxplot(data["after_line"], positions=data['X']+2*wd+0.5, widths=wd)
ax1.boxplot(data["after_alu"], positions=data['X']+3*wd+0.5, widths=wd)
#ax2.boxplot(data["test_coverage"], positions=data['X']+1.75*wd, widths=0.35, patch_artist=True,
#            capprops=dict(color='b'),
#            flierprops=dict(markeredgecolor='b'),
#            medianprops=dict(color='b'))


ax1.set_xticks(data['X']+1.75*wd)
ax1.set_xticklabels(data['XString'], fontsize=15)
#plt.set_yticks(fontsize=15)
ax1.set_title('RTD-Control and RTD-Diploid LINE and Alu False Positive Rate', fontsize=20)
ax1.set_xlabel('Read Splits (Average Read Length)', fontsize=17)
ax1.set_ylabel('Insertion Rate per Billion Bases', fontsize=17)
#ax2.set_ylabel('Read Length (bp)', fontsize=17, color='b')
ax1.legend(loc='upper right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Line_alu_fpr.png")
plt.clf()


