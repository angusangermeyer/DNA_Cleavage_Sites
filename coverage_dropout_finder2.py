#!/usr/bin/env python3

"""
A pipeline to compare sequence mapping coverage ratios between control and experimental conditions.
Allows for the identification of site-specific nuclease cleavage locations.
"""

#Imported modules
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import pathlib
from subprocess import call
import os, sys
import statistics
import math
import numpy as np
from scipy.signal import find_peaks


#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2023 Angus Angermeyer"
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Angus Angermeyer"
__email__ = "angus.angermeyer@gmail.com"
__status__ = "Prototype"



#Variable inputs
return_count = 5
peak_width = 50
plot_width = 500
attP = [3985, 4144] #for plasmid ref

# attP = [10244, 10403] #for PLE dORF10 ref

#Static inputs
fastq_path = "fastq_files/"
CWD = str(pathlib.Path.cwd())

control_seq_replicates = [
	(fastq_path+"MN19E_S222_R1_001.fastq.gz",fastq_path+"MN19E_S222_R2_001.fastq.gz"),
	(fastq_path+"MN20E_S223_R1_001.fastq.gz",fastq_path+"MN20E_S223_R2_001.fastq.gz"),
	(fastq_path+"MN21E_S224_R1_001.fastq.gz",fastq_path+"MN21E_S224_R2_001.fastq.gz")]

experimental_seq_replicates = [
	(fastq_path+"MN4E_S44_R1_001.fastq.gz",fastq_path+"MN4E_S44_R2_001.fastq.gz"),
	(fastq_path+"MN8E_S42_R1_001.fastq.gz",fastq_path+"MN8E_S42_R2_001.fastq.gz"),
	(fastq_path+"MN12E_S46_R1_001.fastq.gz",fastq_path+"MN12E_S46_R2_001.fastq.gz")]

PLE_exp = [control_seq_replicates, experimental_seq_replicates]
exp_names = ["Control", "Experiment"]
fna_path = "pKL06.2_fasta_from_AA_attP.fasta"



fna_record = SeqIO.read(fna_path, "fasta")





call(f"bowtie2-build {fna_path} temp_index 1>/dev/null", shell = True)

name_list = []

for experiment, exp_name in zip(PLE_exp, exp_names):
	raw_values = []
	for replicate in experiment:
	
		forward_fastq = replicate[0]
		reverse_fastq = replicate[1]

		name = forward_fastq[forward_fastq.rindex("/")+1:forward_fastq.index(".fastq.gz")]
		name_list.append(name)

		call(f"bowtie2 -x temp_index -1 {forward_fastq} -2 {reverse_fastq} -p 4 --end-to-end --very-sensitive --no-unal -S temp.sam", shell =True)

		call(f"samtools view -b temp.sam > temp.bam", shell =True) #all reads
		os.remove("temp.sam")
		
		call(f"samtools sort temp.bam > {name}_main_output.bam", shell =True)

		call(f"samtools index {name}_main_output.bam", shell =True)

		call(f"breseq BAM2COV -b {name}_main_output.bam -f {fna_path} -o {CWD}/{name}_{fna_record.id}_output -p 0 -r {fna_record.id}:1-{len(fna_record.seq)} -t", shell =True)
		os.remove("temp.bam")



##Average replicates
		with open(f"{name}_{fna_record.id}_output.tab", "r") as input:
			input.readline()
			x=[]
			y=[]
			read_count = 0
			for line in input:
				line=line.split()
				x.append(int(line[0]))
				reads = int(line[2])+ int(line[3])
				y.append(reads)
				read_count += reads
			
			
			y[:] = [entry / read_count*100 for entry in y] #Percent total reads
# 			y[:] = [entry / read_count for entry in y] #Raw read counts
			raw_values.append(y)

	y_avg = []
	y_dev = []
	with open(f"{exp_name}_avg_data.tab", "w") as output:
		x=1
		for (a, b, c) in zip(raw_values[0], raw_values[1], raw_values[2]):
			avg = sum([a,b,c])/3
			stdev = statistics.stdev([a,b,c])
			output.write(f"{x}\t{avg}\t{stdev}\n")
			y_avg.append(avg)
			y_dev.append(stdev)
			x+=1



###Ratio calculation
with open(f"{exp_names[0]}_avg_data.tab", "r") as pos_input, open(f"{exp_names[1]}_avg_data.tab", "r") as neg_input:
# with open("PLE_pos_avg_data.tab", "r") as pos_input, open("PLE_neg_avg_data.tab", "r") as neg_input:
	log_ratio_list = []
	inverse_log_ratio_list = []
	for pos_line, neg_line in zip(pos_input, neg_input):
		pos_line = pos_line.split()
		neg_line = neg_line.split()
		
		pos_cov = float(pos_line[1])
		neg_cov = float(neg_line[1])
		
		
		#When using an artificial reference with poly-A region where no mapping occurs, the ratio thing fails.
		#Definitely need a better fix, but for the moment I'll just ignore those reasons.
		try:
			ratio = neg_cov/pos_cov
			log_ratio = math.log(ratio,2)
		except ZeroDivisionError:
			log_ratio = 0
			
		log_ratio_list.append(log_ratio)


		try:
			inverse_ratio = pos_cov/neg_cov
			inverse_log_ratio = math.log(inverse_ratio,2)
		
		except (ZeroDivisionError, ValueError) as error:
# 		except ZeroDivisionError:
			inverse_log_ratio = 0

		inverse_log_ratio_list.append(inverse_log_ratio)



with open(f"{exp_names[0]}-{exp_names[1]}_Hits={return_count}_PeakWidth={peak_width}_log2_ratio.txt", "w") as output:
	x=1
	for entry in log_ratio_list:
		output.write(f"{x}\t{entry}\n")
		x+=1



###Peak finding method
peaks = find_peaks(inverse_log_ratio_list, 
			height=None,
			threshold=None,
			distance=peak_width/2,
			prominence=.2,
			width=peak_width,
			wlen=None,
			rel_height=0.5,
			plateau_size=None)[0]

window_list = [[1,1] for x in range(return_count)]
for location in peaks:
	valley = log_ratio_list[location]

	pos = 0
	for entry in window_list:
		if valley < entry[1]:
			window_list.insert(pos,[location, valley])
			window_list = window_list[:-1]
			break
		pos+=1

#Only to force the correct site at t=8
# window_list.insert(0,[10328, -3])

print(window_list)



###Steepest basepair change method
# window_strt = 0
# window_stop = window_strt + 3
# window_list = [[1,0] for x in range(return_count)]
# while window_stop < len(log_ratio_list):
# 	bp_delta = abs(log_ratio_list[window_strt]-log_ratio_list[window_stop])
# 	pos = 0
# 	for item in window_list:
# 		if bp_delta>item[1]:
# 			window_list.insert(pos,[window_strt,bp_delta])
# 			window_list = window_list[:-1]
# 			break
# 		pos+=1
# 		
# 	window_strt += 1
# 	window_stop += 1
# 
# print(window_list)



#Sliding window method. Doesn't work due to overlapping valley windows. Could be fixed, but others seem better.
# window_strt = 0
# window_stop = window_strt + window_size
# window_list = [[1,1] for x in range(return_count)]
# while window_stop <= len(log_ratio_list):
# 	average = sum(log_ratio_list[window_strt:window_stop])/len(log_ratio_list[window_strt:window_stop])
# 	pos = 0
# 	for item in window_list:
# 		if average<item[1]:
# 			window_list.insert(pos,[window_strt,average])
# 			window_list = window_list[:-1]
# 			break
# 		pos+=1
# 		
# 	window_strt += 1
# 	window_stop += 1
# 
# print(window_list)




###Plotting


fig, axs = plt.subplots(len(window_list), figsize=(10, 10))
# fig.suptitle('Vertically stacked subplots')

axis_num = 0
for entry in window_list:
	if entry == [1,1] or len(log_ratio_list)-entry[0]<plot_width/2:

		loc = entry[0]
		x = [i + (1-int(plot_width/2)) for i in range(plot_width)]
		y = [np.nan for i in range(plot_width)]

		axs[axis_num].text(1, -1, "No Peak Found", horizontalalignment='center', verticalalignment='center',)
		axs[axis_num].get_xaxis().set_visible(False)
		
		
		
	else:
		loc = entry[0]
		x = [i + (loc-int(plot_width/2)) for i in range(plot_width)]

		y = []
		for nuc_pos in x:
			try:
				y.append(log_ratio_list[nuc_pos])
			except IndexError:
				y.append(np.nan)
		
		axs[axis_num].plot(x, [0 for i in x], color="grey", linewidth=0.5)
		
		
		if x[0] <= attP[0] and x[-1]>= attP[1]:
			axs[axis_num].axvspan(attP[0], attP[1], ymin=0, ymax=1, facecolor='black', alpha=0.25)
		elif x[-1] > attP[0] and x[-1] < attP[1]:
			axs[axis_num].axvspan(attP[0], x[-1], ymin=0, ymax=1, facecolor='black', alpha=0.25)
		elif x[0] < attP[1] and x[0] > attP[0]:
			axs[axis_num].axvspan(x[0], attP[1], ymin=0, ymax=1, facecolor='black', alpha=0.25)
		else:
			pass


	axs[axis_num].plot(x, y, color="black", linewidth=2)
	axs[axis_num].set_xlim([x[0],x[-1]])
	axs[axis_num].set_ylim(-3,1)

	

	axis_num+=1



fig.set_size_inches(8.5, 2*len(window_list))
plt.tight_layout(pad=2, w_pad=1, h_pad=-1.5)
matplotlib.pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.2)
plt.savefig(f"{exp_names[0]}-{exp_names[1]}_Hits={return_count}_PeakWidth={peak_width}.png", dpi = 600)
plt.close()


