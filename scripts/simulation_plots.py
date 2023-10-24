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
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

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
args = parser.parse_args()

data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break

i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X1.5"].append(float(i)+1.5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 2.5


fig, ax1 = plt.subplots()
wd = 0.5
ax1.bar(data['X'],data["somrit_mean"], color='c', width=wd, edgecolor='k', label="Somrit", alpha=0.4)
ax1.bar(data['X1'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
ax1.bar(data['X2'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2", alpha=0.4)
#ax1.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
ax1.boxplot(data["somrit"], positions=data['X'], widths=wd)
ax1.boxplot(data["tldr"], positions=data['X1'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
#ax1.boxplot(data["xtea"], positions=data['X3'], widths=wd)
ax1.boxplot(data["sniffles"], positions=data['X2'], widths=wd)
ax1.legend(loc='lower right', fontsize=15)

ax1.set_xticks(data['X1'])
ax1.set_xticklabels(data['XString'], fontsize=15)
ax1.tick_params(axis='both', which='minor', size=15)
ax1.set_title('Somrit, tldr and Sniffles2 Recall on Simulated Insertions', fontsize=20)
ax1.set_xlabel('Sample', fontsize=17)
ax1.set_ylabel('Recall', fontsize=17)
ax1.set_ylim([0,1.05])

ax2 = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.75,0.75,0.2,0.2])
ax2.set_axes_locator(ip)
ax2.bar(data['X'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
ax2.boxplot(data["somrit_read"], positions=data['X'], widths=wd)
ax2.set_xticks(data['X'])
ax2.set_xticklabels(data['XString'])
ax2.legend(loc=0)
ax2.set_ylim([0,1.05])

fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_recall.png")
plt.clf()

#exit(0)

data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_500_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_500_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break


i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X4"].append(float(i)+4*wd)
    data["X5"].append(float(i)+5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 3

wd = 0.5
plt.bar(data['X'],data["somrit_mean"], color='r', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.bar(data['X1'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
plt.bar(data['X2'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr Recall", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
plt.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
plt.bar(data['X4'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2 Recall", alpha=0.4)
plt.boxplot(data["somrit"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_read"], positions=data['X1'], widths=wd)
plt.boxplot(data["tldr"], positions=data['X2'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
plt.boxplot(data["xtea"], positions=data['X3'], widths=wd)
plt.boxplot(data["sniffles"], positions=data['X4'], widths=wd)

plt.xticks(data['X2'], data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, tldr, xTea-Long and Sniffles2 Recall on Simulated Insertions', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='lower right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_500_recall.png")
plt.clf()

data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_2000_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_2000_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break


i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X4"].append(float(i)+4*wd)
    data["X5"].append(float(i)+5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 3

wd = 0.5
plt.bar(data['X'],data["somrit_mean"], color='r', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.bar(data['X1'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
plt.bar(data['X2'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr Recall", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
plt.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
plt.bar(data['X4'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2 Recall", alpha=0.4)
plt.boxplot(data["somrit"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_read"], positions=data['X1'], widths=wd)
plt.boxplot(data["tldr"], positions=data['X2'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
plt.boxplot(data["xtea"], positions=data['X3'], widths=wd)
plt.boxplot(data["sniffles"], positions=data['X4'], widths=wd)

plt.xticks(data['X2'], data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, tldr, xTea-Long and Sniffles2 Recall on Simulated Insertions', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='lower right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_2000_recall.png")
plt.clf()


data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_6000_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_6000_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break


i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X4"].append(float(i)+4*wd)
    data["X5"].append(float(i)+5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 3

wd = 0.5
plt.bar(data['X'],data["somrit_mean"], color='r', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.bar(data['X1'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
plt.bar(data['X2'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr Recall", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
plt.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
plt.bar(data['X4'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2 Recall", alpha=0.4)
plt.boxplot(data["somrit"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_read"], positions=data['X1'], widths=wd)
plt.boxplot(data["tldr"], positions=data['X2'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
plt.boxplot(data["xtea"], positions=data['X3'], widths=wd)
plt.boxplot(data["sniffles"], positions=data['X4'], widths=wd)

plt.xticks(data['X2'], data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, tldr, xTea-Long and Sniffles2 Recall on Simulated Insertions', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='lower right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_6000_recall.png")
plt.clf()


data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_LINE_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_LINE_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break


i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X4"].append(float(i)+4*wd)
    data["X5"].append(float(i)+5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 3

wd = 0.5
plt.bar(data['X'],data["somrit_mean"], color='r', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.bar(data['X1'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
plt.bar(data['X2'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr Recall", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
plt.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
plt.bar(data['X4'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2 Recall", alpha=0.4)
plt.boxplot(data["somrit"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_read"], positions=data['X1'], widths=wd)
plt.boxplot(data["tldr"], positions=data['X2'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
plt.boxplot(data["xtea"], positions=data['X3'], widths=wd)
plt.boxplot(data["sniffles"], positions=data['X4'], widths=wd)

plt.xticks(data['X2'], data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, tldr, xTea-Long and Sniffles2 Recall on Simulated Insertions', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='lower right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_LINE_recall.png")
plt.clf()




data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_Alu_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_Alu_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break


i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X4"].append(float(i)+4*wd)
    data["X5"].append(float(i)+5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 3

wd = 0.5
plt.bar(data['X'],data["somrit_mean"], color='r', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.bar(data['X1'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
plt.bar(data['X2'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr Recall", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
plt.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
plt.bar(data['X4'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2 Recall", alpha=0.4)
plt.boxplot(data["somrit"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_read"], positions=data['X1'], widths=wd)
plt.boxplot(data["tldr"], positions=data['X2'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
plt.boxplot(data["xtea"], positions=data['X3'], widths=wd)
plt.boxplot(data["sniffles"], positions=data['X4'], widths=wd)

plt.xticks(data['X2'], data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, tldr, xTea-Long and Sniffles2 Recall on Simulated Insertions', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='lower right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Alu_recall.png")
plt.clf()



data = defaultdict(list)
# Read and Plot Sniffles/CuteSV Before after Precision & Recall
xtea_recall = defaultdict(list)
sniffles_recall = defaultdict(list)
tldr_recall = defaultdict(list)
somrit_recall = defaultdict(list)
tldr_read_recall = defaultdict(list)
somrit_read_recall = defaultdict(list)
for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_results_SVA_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_results_SVA_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_tsv:
            for line in in_tsv:
                row = line.strip().split("\t")
                somrit_recall[sample].append(float(row[3]))
                somrit_read_recall[sample].append(float(row[6]))
                tldr_recall[sample].append(float(row[9]))
                tldr_read_recall[sample].append(float(row[12]))
                xtea_recall[sample].append(float(row[15]))
                sniffles_recall[sample].append(float(row[18]))
                break


i = 1
wd = 0.5
count = 0
for sample in args.samples.split(','):
    data["X"].append(float(i))
    data["X1"].append(float(i)+wd)
    data["X2"].append(float(i)+2*wd)
    data["X3"].append(float(i)+3*wd)
    data["X4"].append(float(i)+4*wd)
    data["X5"].append(float(i)+5*wd)
    data["XString"].append(sample)
    data["sniffles_mean"].append(np.mean(sniffles_recall[sample]))
    data["sniffles"].append(sniffles_recall[sample])
    data["somrit_mean"].append(np.mean(somrit_recall[sample]))
    data["somrit"].append(somrit_recall[sample])
    data["somrit_read_mean"].append(np.mean(somrit_read_recall[sample]))
    data["somrit_read"].append(somrit_read_recall[sample])
    data["xtea_mean"].append(np.mean(xtea_recall[sample]))
    data["xtea"].append(xtea_recall[sample])
    data["tldr_mean"].append(np.mean(tldr_recall[sample]))
    data["tldr"].append(tldr_recall[sample])
    data["tldr_read_mean"].append(np.mean(tldr_read_recall[sample]))
    data["tldr_read"].append(tldr_read_recall[sample])
    i += 3

wd = 0.5
plt.bar(data['X'],data["somrit_mean"], color='r', width=wd, edgecolor='k', label="Somrit Recall", alpha=0.4)
plt.bar(data['X1'],data["somrit_read_mean"], color='y', width=wd, edgecolor='k', label="Somrit Read Recall", alpha=0.4)
plt.bar(data['X2'],data["tldr_mean"], color='b', width=wd, edgecolor='k', label="tldr Recall", alpha=0.4)
#plt.bar(data['X3'],data["tldr_read_mean"], color='g', width=wd, edgecolor='k', label="tldr Read Recall", alpha=0.4)
plt.bar(data['X3'],data["xtea_mean"], color='c', width=wd, edgecolor='k', label="xTea-Long Recall", alpha=0.4)
plt.bar(data['X4'],data["sniffles_mean"], color='m', width=wd, edgecolor='k', label="Sniffles2 Recall", alpha=0.4)
plt.boxplot(data["somrit"], positions=data['X'], widths=wd)
plt.boxplot(data["somrit_read"], positions=data['X1'], widths=wd)
plt.boxplot(data["tldr"], positions=data['X2'], widths=wd)
#plt.boxplot(data["tldr_read"], positions=data['X3'], widths=wd)
plt.boxplot(data["xtea"], positions=data['X3'], widths=wd)
plt.boxplot(data["sniffles"], positions=data['X4'], widths=wd)

plt.xticks(data['X2'], data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Somrit, tldr, xTea-Long and Sniffles2 Recall on Simulated Insertions', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Recall', fontsize=17)
plt.legend(loc='lower right', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_SVA_recall.png")
plt.clf()






sniffles_before_mem = defaultdict(list)
sniffles_before_time = defaultdict(list)
tldr_before_mem = defaultdict(list) 
tldr_before_time = defaultdict(list)
xtea_before_mem = defaultdict(list)
xtea_before_time = defaultdict(list)
somrit_mem = defaultdict(list)
somrit_time = defaultdict(list)
somrit_ideal_mem = defaultdict(list)
somrit_ideal_time = defaultdict(list)
somrit_realign_only_mem = defaultdict(list)
somrit_realign_only_time = defaultdict(list)


for sample in args.samples.split(','):
    for rep in [1,2,3,4,5,6,7,8,9,10,11,12]:
        # check if the file exists
        if not os.path.isfile(args.input_dir+"/"+sample+"/simulation_benchmark_"+sample.replace('0','')+".spike_in."+str(rep)+".txt"):
            continue
        with open(args.input_dir+"/"+sample+"/simulation_benchmark_"+sample.replace('0','')+".spike_in."+str(rep)+".txt", 'r') as in_time_mem:
            for line in in_time_mem:
                row = line.strip().split("\t")
                sniffles_before_mem[sample].append(float(row[5]))
                sniffles_before_time[sample].append(float(row[6]))
                xtea_before_mem[sample].append(float(row[7]))
                xtea_before_time[sample].append(float(row[8]))
                tldr_before_mem[sample].append(float(row[9]))
                tldr_before_time[sample].append(float(row[10]))
                somrit_extract_mem = float(row[11])
                somrit_extract_time= float(row[12])
                somrit_realign_mem = float(row[13])
                somrit_realign_time= float(row[14])
                somrit_realign_max_mem = float(row[15])
                somrit_realign_max_time= float(row[16])
                somrit_realign_avg_mem = float(row[17])
                somrit_realign_avg_time= float(row[18])
                somrit_merge_mem = float(row[19])
                somrit_merge_time= float(row[20])
                somrit_classify_mem = float(row[21])
                somrit_classify_time= float(row[22])
                somrit_filter_mem = float(row[23])
                somrit_filter_time= float(row[24])
                somrit_time[sample].append(somrit_extract_time + somrit_realign_time + somrit_merge_time + somrit_classify_time + somrit_filter_time)
                somrit_mem[sample].append(max(somrit_extract_mem, somrit_realign_max_mem, somrit_merge_mem, somrit_classify_mem, somrit_filter_mem))
                somrit_realign_only_mem[sample].append(somrit_realign_max_mem)
                somrit_realign_only_time[sample].append(somrit_realign_time)
                somrit_ideal_mem[sample].append(max(somrit_extract_mem, somrit_realign_max_mem, somrit_merge_mem, somrit_classify_mem, somrit_filter_mem))
                somrit_ideal_time[sample].append(somrit_extract_time + somrit_realign_max_time + somrit_merge_time + somrit_classify_time + somrit_filter_time)

print("data")

# Compute the mean and stdev of each measure
data = defaultdict(list)
count = 0
for sample in args.samples.split(','):
    data['X'].append(count*3)
    data['XString'].append(sample)
    data["sniffles_before_mem_mean"].append(np.mean(sniffles_before_mem[sample]))
    data["sniffles_before_time_mean"].append(np.mean(sniffles_before_time[sample]))
    data["xtea_before_mem_mean"].append(np.mean(xtea_before_mem[sample]))
    data["xtea_before_time_mean"].append(np.mean(xtea_before_time[sample]))
    data["tldr_before_mem_mean"].append(np.mean(tldr_before_mem[sample]))
    data["tldr_before_time_mean"].append(np.mean(tldr_before_time[sample]))
    data["somrit_mem_mean"].append(np.mean(somrit_mem[sample]))
    data["somrit_time_mean"].append(np.mean(somrit_time[sample]))
    data["somrit_realign_only_mem_mean"].append(np.mean(somrit_realign_only_mem[sample]))
    data["somrit_realign_only_time_mean"].append(np.mean(somrit_realign_only_time[sample]))
    data["somrit_ideal_mem_mean"].append(np.mean(somrit_ideal_mem[sample]))
    data["somrit_ideal_time_mean"].append(np.mean(somrit_ideal_time[sample]))
    data["sniffles_before_mem"].append(sniffles_before_mem[sample])
    data["sniffles_before_time"].append(sniffles_before_time[sample])
    data["xtea_before_mem"].append(xtea_before_mem[sample])
    data["xtea_before_time"].append(xtea_before_time[sample])
    data["tldr_before_mem"].append(tldr_before_mem[sample])
    data["tldr_before_time"].append(tldr_before_time[sample])
    data["somrit_mem"].append(somrit_mem[sample])
    data["somrit_time"].append(somrit_time[sample])
    data["somrit_ideal_mem"].append(somrit_ideal_mem[sample])
    data["somrit_ideal_time"].append(somrit_ideal_time[sample])
    data["somrit_realign_only_mem"].append(somrit_realign_only_mem[sample])
    data["somrit_realign_only_time"].append(somrit_realign_only_time[sample])
    count += 1

data['X'] = np.array(data['X'])
# Plot Memory Graph
wd = 0.4
plt.bar(data['X']+wd,data["sniffles_before_mem_mean"], color='m', alpha=0.4, width=wd, edgecolor='k', label="Sniffles")
plt.bar(data['X']+2*wd,data["xtea_before_mem_mean"], color='r', alpha=0.4, width=wd, edgecolor='k', label="xTea-Long")
plt.bar(data['X']+3*wd,data["tldr_before_mem_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="tldr")
plt.bar(data['X']+4*wd,data["somrit_mem_mean"], color='c', alpha=0.4, width=wd, edgecolor='k', label="somrit")
plt.boxplot(data["sniffles_before_mem"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["xtea_before_mem"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["tldr_before_mem"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_mem"], positions=data['X']+4*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Max Memory Usage per Tool', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Memory Usage (MiB)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Memory.png")
plt.clf()


plt.bar(data['X']+wd,data["sniffles_before_time_mean"], color='m', alpha=0.4, width=wd, edgecolor='k', label="Sniffles")
plt.bar(data['X']+2*wd,data["xtea_before_time_mean"], color='r', alpha=0.4, width=wd, edgecolor='k', label="xTea-Long")
plt.bar(data['X']+3*wd,data["tldr_before_time_mean"], color='b', alpha=0.4, width=wd, edgecolor='k', label="tldr")
plt.bar(data['X']+4*wd,data["somrit_ideal_time_mean"], alpha=0.4, color='c', width=wd, edgecolor='k', label="somrit")
plt.boxplot(data["sniffles_before_time"], positions=data['X']+wd, widths=wd)
plt.boxplot(data["xtea_before_time"], positions=data['X']+2*wd, widths=wd)
plt.boxplot(data["tldr_before_time"], positions=data['X']+3*wd, widths=wd)
plt.boxplot(data["somrit_ideal_time"], positions=data['X']+4*wd, widths=wd)

plt.xticks(data['X']+2.5*wd, data['XString'], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Max Time per Tool', fontsize=20)
plt.xlabel('Sample', fontsize=17)
plt.ylabel('Time (s)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_Time.png")
plt.clf()

