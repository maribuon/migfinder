# 2016.01.20
# Metagenomic Integron-associated Gene finder

# 1) uses HattCI to predict attC sites
# 2) uses Infernal to validate attC sites secondary structure
# 3) uses Prodigal to annotate CDS associated with the attC-sites found


# Requires:
	# HattCI, Infernal or Prodigal to be installed

# Filtering results, i.e. only what has been picked by HattCI and then scored by Infernal
# Infernal @ max


#---------------------------------------------------------------------------#
import os
import re
import csv
import sys
import argparse
import subprocess

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

from Bio import SeqIO

import logging

#---------------------------------------------------------------------------#
# HattCI call
def hattci(fastafile, output_directory, both=True, nseq=1000, nthread=6):
	# creating a dir for the hmm results
	hmmer_results = os.path.join(output_directory, "hmmresults")
	os.makedirs(hmmer_results, exist_ok=True)
	
	# os.chdir("hmmresults")
	output_file = f"{hmmer_results}/{os.path.basename(fastafile)}"
	output_file_tmp = f"{output_file}.tmp"
	output_file_log = f"{output_file}.hmmlog"
	
	# calling hattci, both strands or not?
	params = [
		"hattci.out",
		"-s",
		str(nseq),
		"-t",
		str(nthread),
		fastafile,
		output_file_tmp
	]
	if both:
		params.insert(1, "-b")
	
	out_f=open(output_file_log, "w")
	subprocess.run(params, stdout=out_f)

	#--------------#
	# parsing file
	#--------------#
	with open(output_file_tmp, "r") as file_in, open(output_file + ".hmm", "w") as file_out:
		# extracting only table from the results
		for line in file_in.readlines():
			# removing first line with column names
			if 'hit' in line:
				continue
			# stopping if ------------- line is found
			elif re.match(r'-----',line):
				break
			else:
				file_out.write(line)
	os.unlink(output_file_tmp)
	# TODO: it is possible that outHattCI.fasta will be written in the base dir
	os.rename(f"{output_directory}/outHattCI.fasta", f"{output_file}_hattci.fasta")
	
	return f"{output_file}_hattci.fasta"

#---------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
# Infernal call
# NOTE: Infernal score both strands and do not pick one.
# HattCI picks the best
# therefore, when results are combined, only one strand can be selected.
def infernal(fastafile, output_directory, cm_model):
	# creating a dir for the hmm results
	cmresults_out_dir = f"{output_directory}/cmresults"
	output_file = os.path.basename(fastafile)
	output_file_tmp = f"{cmresults_out_dir}/{output_file}.tmp"

	os.makedirs(cmresults_out_dir, exist_ok=True)
	params = [
		"cmsearch",
		"--max",
		"-o",
		output_file_tmp,
		cm_model,
		fastafile
	]
	# calling infernal
	subprocess.run(params)

	#--------------#
	# parsing file
	#--------------#
	with open(output_file_tmp) as filein, open(output_file + ".cm",'w') as fileout, open(output_file + ".cmlog",'w') as filelog:
		# extracting only table from the results
		for line in filein.readlines():
			# creating cmlog
			if re.match(r'\#',line):
				filelog.write(line)	
			# writing the table part
			elif re.search('\(\d+\)',line):
				fileout.write(line)
			elif "Hit alignments" in line:
				break
	os.unlink(output_file_tmp)
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Filtering: Process CM results into only one table file
def posproc(fastafile, output_directory, k_cm=20, dist_threshold=4000):
	remove_attC_wrongstrand = 0
	# --------------- CM -------------- #
	# converting cm into an array
	# extracting tags: acc number + sta + sto
	# NOTE: because we are filtering HattCI -> Infernal, infernal tag is acc_sta_sto, have to split("_") to extract tags
	# extracting score	--> cmscore
	# extracting e-value	--> cmevalue
	tmp = "cmresults/"+out+".cm"
	fcm = list(csv.reader(open(tmp, 'rb'),delimiter='\t'))
	Mcm = len(fcm)
	cmtag = []
	cmscore = []
	cmevalue = []
	sta = []
	sto = []
	names = []
	fhmm = []		# Viterbi score from hmm (Vscore)
	table = []
	strand = []
	tag = []
	m = 0
	while m < Mcm:
		fcm[m] = fcm[m][0].split()
		# removing hits with score lower than threshold k_cm
		# since CM results are sorted, breaking is sufficient to remove all <k_cm	
		if float(fcm[m][3]) < k_cm:
			break			
		# remove end
		else:
			tag.append(fcm[m][5])
			aux = fcm[m][5].split("_")
			names.append("_".join(aux[0:len(aux)-2]))
			cmtag.append(names[m]+"_"+aux[len(aux)-2]+'_'+aux[len(aux)-1])
			aux_sta = int(aux[-2])
			aux_sto = int(aux[-1])
			cmscore.append(fcm[m][3])
			cmevalue.append(fcm[m][2])
		# The section below will use coordinates from HattCI, which is desirable to maintain R'' and R' sites.
			if aux_sta <= aux_sto:
				sta.append(aux_sta)
				sto.append(aux_sto)
				strand.append("+")
			else:
				sto.append(aux_sta)
				sta.append(aux_sto)
				strand.append("-")
			last = fcm[m][len(fcm[m])-1]
			if last == 'REVERSED':
				fhmm.append(fcm[m][len(fcm[m])-2])
			else:
				fhmm.append(fcm[m][len(fcm[m])-1])
			# saving table
			table.append([names[m],int(sta[m]),sto[m], strand[m],fhmm[m],cmscore[m],cmevalue[m]])		
		m = m + 1
	# ----------------- sorting ----------------------- #
	table = sorted(table, key = lambda data:data[1])
	table = sorted(table, key = lambda data:data[0])
	# ------------------------------------------------- #
	# Infernal will evaluate and maybe pick both strands: Redundancy
	# removing Redundancy and attC site on the wrong strand in the same an integron
	# 	we are selecting HattCI coord, so it may happen that we have double entry,
	#	i.e. if CM found the same hit on top and bottom strand, remove one with lower CMscore
	Mcm = len(table)
	#if Mcm > 0
	m = 1
	while m < Mcm:
		if table[m][0] == table[m-1][0] and table[m][1] == table[m-1][1]:
			if float(table[m][5]) >= float(table[m-1][5]):
				del table[m-1]
			else:
				del table[m]
			m -= 1
		m += 1
		Mcm = len(table)
	# removing wrong strand and other redundancy
	m = 0
	while m < Mcm-1:
		# counting no. attC sites on each strand
		n_plus = 0
		n_minus = 0
		# saving the indeces for minus and plus
		ind_plus = []
		ind_minus = []
		# saving scores, as cumsum
		score_plus = 0
		score_minus = 0
		if table[m][3] == '+':
			n_plus = n_plus + 1
			ind_plus.append(m)
			score_plus = score_plus + float(table[m][4])
		else:
			n_minus = n_minus + 1
			ind_minus.append(m)
			score_minus = score_minus + float(table[m][4])
		# end counting
		# --- picking hits on the same integron
		n = m + 1
		while table[m][0] == table[n][0] and int(table[n][1])-int(table[n-1][2]) <= dist_threshold:
			# start counting strands
			if table[n][3] == '+':
				n_plus = n_plus + 1
				ind_plus.append(n)
				score_plus = score_plus + float(table[n][5])
			else:
				n_minus = n_minus + 1
				ind_minus.append(n)
				score_minus = score_minus + float(table[n][5])
			# end counting strands
			n = n + 1
			if n >= Mcm:
				break
		# --- end picking same integron
		# ---- start removing "wrong" strand:
		if n_plus > n_minus:
			if n_minus > 0:
				for i in range(n_minus-1,-1,-1):
					del table[ind_minus[i]]
					remove_attC_wrongstrand += 1
					Mcm = Mcm - 1
				# end for
			n_remove = n_minus
		elif n_minus > n_plus:
			if n_plus > 0:
				for i in range(n_plus-1,-1,-1):
					del table[ind_plus[i]]
					remove_attC_wrongstrand += 1
					Mcm = Mcm - 1
				# end for
			n_remove = n_plus
		elif n_plus == n_minus:
			if score_plus >= score_minus:
				for i in range(n_minus-1,-1,-1):
					del table[ind_minus[i]]
					remove_attC_wrongstrand += 1
					Mcm = Mcm - 1
				# end for
			else:
				for i in range(n_plus-1,-1,-1):
					del table[ind_plus[i]]
					remove_attC_wrongstrand += 1
					Mcm = Mcm - 1
				# end for
			n_remove = n_minus
		# --- end removing "wrong" strand
		m = n - n_remove
		Mcm = len(table)
	# ----------------- Removing singletons ----------------------- #
	Mcm = len(table)
	m = 0
	# dealing with the first position, just check coming one.
	if len(table) > 1:
		while m == 0:
			if m+1 < Mcm and table[m][0] != table[m+1][0]:
				del table[m]
				m -= 1
			elif m+1 < Mcm and table[m][0] == table[m+1][0] and table[m+1][1] - table[m][2] > dist_threshold:
				del table[m]
				m -= 1
			m += 1
			Mcm = len(table)
	else:
		table = []
	# dealing with the other positions, check previous and next
	if len(table) > 1:
		m = 1
		while m < Mcm-1:
			# singleton
			if table[m][0] != table[m-1][0] and table[m][0] != table[m+1][0]:
				del table[m]
				m -= 1
			# singleton after an integron (same seq) [(no other integron) or (a third integron)]
			elif (table[m][0] == table[m-1][0] and table[m][0] != table[m+1][0] and  table[m][1] - table[m-1][2] > dist_threshold) or (table[m][0] == table[m-1][0] and table[m][0] == table[m+1][0] and  table[m][1] - table[m-1][2] > dist_threshold and table[m+1][1] - table[m][2] > dist_threshold):
				del table[m]
				m -= 1
			# singleton before an integron (same seq) [(no other integron) or (a third integron)]
			elif (table[m][0] != table[m-1][0] and table[m][0] == table[m+1][0] and  table[m+1][1] - table[m][2] > dist_threshold) or (table[m][0] == table[m-1][0] and table[m][0] == table[m+1][0] and table[m+1][1] - table[m][2] > dist_threshold and table[m][1] - table[m-1][2] > dist_threshold):
				del table[m]
				m -= 1
			m += 1
			Mcm = len(table)
		if Mcm > 1 and table[Mcm-1][0] != table[Mcm-2][0]:
			del table[Mcm-1]
		elif Mcm > 1 and table[Mcm-1][0] == table[Mcm-2][0] and table[Mcm-1][1] - table[Mcm-2][2] > dist_threshold:
			del table[Mcm-1]
		Mcm = len(table)
	else:
		table = []
	#
	logging.info(f"attC hits from pipeline: {len(table)}")
	with open(out+".filtering",'w') as fileout:
		fileout.write("Total_attC:	" + str(len(table)) + '\n')
		fileout.write("Del_attC_ws:	" + str(remove_attC_wrongstrand) + '\n')
	# ----------- Saving results table -------------- #
	Mcm = len(table)
	if table:
		# selecting coordinates for the table from hmmresults to cmresults
		tmp = "hmmresults/" + out + ".hmm"
		fhmm = list(csv.reader(open(tmp, 'rb'),delimiter='\t'))
		Mhmm = len(fhmm)
		# opening file to save
		fileout = open("cmresults/"+out+"_attC.res",'w')
		names_sel = []
		for m in range(0,Mcm):
			names_sel.append(table[m][0])
			table[m][1] = str(table[m][1])
			table[m][2] = str(table[m][2])
			for n in range(0,Mhmm):
				# top strand
				if table[m][0] == fhmm[n][1].strip() and table[m][1] == fhmm[n][2].strip() and table[m][2] == fhmm[n][3].strip():
					kk = 7
					for k in range(5,12):
						table[m].insert(kk,fhmm[n][k].strip())
						kk += 1
					break
				# bottom strand
				elif table[m][0] == fhmm[n][1].strip() and table[m][1] == fhmm[n][3].strip() and table[m][2] == fhmm[n][2].strip():
					kk = 7
					for k in range(5,12):
						table[m].insert(kk,fhmm[n][k].strip())
						kk += 1
					break
			fileout.write("\t".join(table[m]))
			fileout.write("\n")
		fileout.close()
	# ------------ Selecting only sequences with attC hit to predict ORF -----#
	if table:
		names_sel = set(names_sel)
		names_sel = sorted(names_sel)
		all_sequences = SeqIO.parse(open(fastafile), "fasta")
		hits_sequences = "cmresults/"+out+"_infernal.fasta"
		with open(hits_sequences, "w") as f:
			for seq in all_sequences:
				if seq.id in names_sel:
					foo = SeqIO.write([seq],f,"fasta")
		f.close()
	#
#---------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
# Prodigal call
# therefore, when results are combined, only one strand can be selected.
def prodigal(fastafile, output_directory, save_orf):
	# creating a dir for the prodigal results
	orfresults_out_dir = f"{output_directory}/orfresults"
	output_file = os.path.basename(fastafile)
	output_file_gff = f"{orfresults_out_dir}/{output_file}_orf.gff"
	
	os.makedirs(orfresults_out_dir, exist_ok=True)
	
	params = [
		"prodigal",
		"-i",
		fastafile,
		"-p",
		"meta",
		"-o",
		output_file_gff,
		"-f",
		"gff",
		"-q",
		"-c"
	]
	
	if save_orf:
		output_file_faa = f"{orfresults_out_dir}/{output_file}_ALLorf.faa"
		output_file_fna = f"{orfresults_out_dir}/{output_file}_ALLorf.fna"	
		additional_param=[
			"-a",
			output_file_faa,
			"-d",
			output_file_fna
		]
		params=params+additional_param

	# calling prodigal
	subprocess.run(params)

	# Parsing the gff output file to extract CDS
	with open(output_file_gff) as gffin, open("tmp.gff",'w') as fileout:
		for line in gffin.readlines():
			if line.split('\t')[2] == 'CDS':
				fileout.write(line)
				
	os.rename("tmp.gff", f"{output_file_gff}")

#---------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
# Filtering2: Process CM results into only one table file
def posproc2(output_directory, k_orf = 0, d_CDS_attC = 500, dist_threshold=4000):
	remove_not_int = 0
	remove_k = 0
	remove_overlap = 0
	remove_overlap_attC = 0
	# -------------- opening attC predictions ------------------- #
	tmp = "cmresults/"+out+"_attC.res"
	attc = list(csv.reader(open(tmp, 'rb'),delimiter='\t'))
	MattC = len(attc)
	data = []
	# editing the attC_result matrix
	for m in range(0,MattC):
		data.append([attc[m][0], "attC_site", int(attc[m][1]), attc[m][2], attc[m][3], attc[m][5], int(attc[m][2])-int(attc[m][1])+1, attc[m][4], attc[m][7], attc[m][8], attc[m][9], attc[m][10], attc[m][11], attc[m][12], attc[m][13]   ])
	data = sorted(data, key = lambda data:data[2])
	data = sorted(data, key = lambda data:data[0])
	m = 1		
	while m < MattC:
		if data[m-1] == data[m]:
			logging.warning("Warning: duplicated hit (attC)! Check for duplicated entry in the fasta file")
			del data[m]
			m -= 1
			MattC = len(data)
		m += 1
	MattC = len(data)
	# ----------------------------------------------------------- #
	# numbering attC_sites on the same integron
	# number attC_#, where # is the number of attC sites on the same integron
	hits = 0
	nintegron = 0
	if MattC > 1:
		m = 0
		while m < (MattC):
			# first is a flag to say if the first attC site (when >1 on the same sequence has already been written to a file)
			first = False
			nattC = 0
			data[m].insert(7,"singleton")
			data[m].insert(8,0)
			tap = data[m][4]
			#
			mm = m +1
			while mm < MattC:
				if data[m][0] == data[mm][0] and data[mm][2]-int(data[mm-1][3])< dist_threshold:
					nintegron = nintegron +1
					if (first == True):
						nattC = nattC+1
						nintegron = nintegron -1
						hits = hits +1
						# write mm
						dist = data[mm][2] - int(data[mm-1][3])
						data[mm].insert(7,"attC_"+str(nattC))
						data[mm].insert(8,dist)
					else:
						nattC = nattC+2
						hits = hits+2
						# write m, re-writting. Singleton removed, attC_1 for the 1st attC on an integron
						data[mm-1][7] = "attC_1"
						data[mm-1][8] = 0
						# write mm
						dist = data[mm][2] - int(data[mm-1][3])
						data[mm].insert(7,"attC_"+str(nattC))
						data[mm].insert(8,dist)
					first = True
					mm = mm +1
				else:
					break
			if nattC > 0 and tap == '-':
				for i in range(0,nattC):
					data[m+i][7] = "attC_"+str(nattC-i)
			m = mm
	# ---------------------------------------------------------- #
	if MattC > 1:
	# -------------- opening ORF predictions ------------------- #
		tmp = f"{orfresults_out_dir}/{output_file}_orf.gff"
		orf = list(csv.reader(open(tmp, 'rb'),delimiter='\t'))
		total_orf = len(orf)
		Morf = len(orf)
		# adding the orf to the data matrix
		for m in range(0,Morf):
			if float(orf[m][5]) > k_orf:
				data.append([orf[m][0], orf[m][2], int(orf[m][3]), orf[m][4], orf[m][6], orf[m][5], int(orf[m][4])-int(orf[m][3])+1, orf[m][8] ])
			else:
				Morf = Morf - 1
				remove_k += 1
		# -------------- sorting ----------------------- #
		data = sorted(data, key = lambda data:data[2])
		data = sorted(data, key = lambda data:data[0])
		# it may happens that there are duplicated sequences in a fasta file, so remove identical hits
		M = len(data)
		m = 1	
		while m < M:
			if data[m-1][0] == data[m][0] and data[m-1][1] == data[m][1] and data[m-1][2] == data[m][2] and data[m-1][3] == data[m][3]:
				logging.warning("Warning: duplicated hit! Check for duplicated entry in the fasta file")
				del data[m]
				m -= 1
				M = len(data)
			m += 1
		M = len(data)
		#	---------------- Removing CDS on the wrong strand --------- #
		# CDS can happen on both strand as long as they do not overlap beyond attC site
		# Therefore they will not be removed if on the wrong strand.
		# This seems to be more common on superintegrons.
		# ------ Calculate dist between elements on the same seq ---- #
		for m in range(0,M):
			# creating distance_between_elements column
			if m >= 1 and data[m][0] == data[m-1][0]:
				data[m].insert(7, data[m][2] - int(data[m-1][3]))
			else:
				data[m].insert(7,"ini0")
		# --------------- Remove overlaping elements --------------#
		# removing distance < -50, else ignore.
		# remove always the CDS.
		m = 0
		while m < M:
			if data[m][7] < -50:
				if data[m][1] == "CDS":
					del data[m]
					remove_overlap += 1
					#m -= 1
					# correcting distances, note now m is updated for a
					if m > 0 and m < len(data) and data[m][0] == data[m-1][0]:
						data[m][7] = data[m][2] - int(data[m-1][3])
					elif m < len(data):
						data[m][7] = "ini0"
					m = m - 1
				elif data[m][1] == "attC_site":
					if data[m-1][1] == "CDS":
						del data[m-1]
						remove_overlap += 1
						# m has to be adjusted before, so that dist can be recalculated
						m = m - 1
						# correcting distances, note now m is updated for a
						if data[m][0] == data[m-1][0]:
							data[m][7] = data[m][2] - int(data[m-1][3])
						else:
							data[m][7] = "ini0"
					elif data[m-1][1] == "attC_site":
						if float(data[m][5]) >= float(data[m-1][5]):
							del data[m-1]
							remove_overlap_attC += 1
							# m has to be adjusted before, so that dist can be rcalculated
							m = m - 1
							# correcting distances, note now m is updated for a
							if data[m][0] == data[m-1][0]:
								data[m][7] = data[m][2] - int(data[m-1][3])
							else:
								data[m][7] = "ini0"
						else:
							del data[m]
							remove_overlap_attC += 1
							# correcting distances, note now m is updated due to deletion
							if m < M - 1:	# if m is the last one no need to correct distances
								if data[m][0] == data[m-1][0]:
									data[m][7] = data[m][2] - int(data[m-1][3])
								else:
									data[m][7] = "ini0"
							m = m - 1
			else:			
				m = m  + 1
			M = len(data)
		# ---- Remove CDSs that are before and beyond attC ---- #
		# Note that: 	On top strand: 		CDS, attC, CDS, attC...
		#		On bottom strand:	attC, CDS, attC, CDS
		m = 0
		while m < (M - 1):
			if data[m][1] == "CDS":
				attCind = []
				# 2+ integrons may happen on the same seq. int_flag keeps track of that
				int_flag = False
				mm = m +1
				while mm < M and data[mm][0] == data[m][0]:
					# if int_flag = True, one int already found
					# if also data[mm][9] == 0, it means a 2nd integron has been found, so break
					if data[mm][1] == "attC_site" and data[mm][9] == 0 and int_flag:
						# mm needs to be updated so that we do not count CDSs on the 2nd int
						mm = attCind[-1] + 1
						break
					if data[mm][1] == "attC_site":
						attCind.append(mm)
						strand = data[mm][4]
						if data[mm][9] == 0:
							int_flag = True
					mm = mm + 1
				# end while
				# if no attC --> remove all, independtly of strand
				if not attCind:
					del data[m:mm]
					remove_not_int += (mm - m)
					mm = mm - (mm-m)
				# if many attC_sites
				else:
					# correcting mm to be the 1st entry after the last attC on the integron.
					last = mm						
					mm = attCind[-1] + 1
					if strand == '+':
						# if more than 1 CDS before the 1st attC_site
						# --> remove all CDS minus 1
						if attCind[0] > (m+1):
							del data[m:(attCind[0]-1)]
							remove_not_int += ((attCind[0]-1) - m)
							mm = mm - ((attCind[0]-1)-m)
							# adjusting dist between elements for the remaining sites
							data[m][7] = "ini0"
							Nremove = (attCind[0]-1)-m
							attCind[:] = [x - Nremove for x in attCind]
							last = last - Nremove
						# if CDSs exist after the last attC_site
						# --> remove all CDS after the last attCsite
						if attCind[-1] < (last-1):
							del data[(attCind[-1]+1):last]
							remove_not_int += (last - (attCind[-1]+1))
							Nremove = (last-(attCind[-1]+1))
							mm = mm - Nremove
							last = last - Nremove
					else:
						remove1 = 0
						remove2 = 0
						# if more than 2+ CDS after the last attC_site
						# --> remove all 2+ CDS
						if attCind[-1] < (last-2):
							del data[(attCind[-1]+2):last]
							remove_not_int += (last - (attCind[-1]+2))
							remove1 = (last-(attCind[-1]+2))
						# if CDSs exist before the 1st attC_site
						# --> remove all CDS before the 1st attCsite
						if attCind[0] > (m):
							del data[m:attCind[0]]
							remove_not_int += (attCind[0]-m)
							remove2 = (attCind[0]-m)
							# adjusting dist between elements for the remaining sites
							data[m][7] = "ini0"
							attCind[:] = [x - (attCind[0]-m) for x in attCind]
						mm = mm - remove1 - remove2
						last = last - remove1 - remove2
				m = last
				M = len(data)
			# if the integron starts by an attC_site
			else:
				attCind = []
				strand = data[m][4]
				# Note: here we are starting by finding an attC_site, thus int_flag is irrelevant!
				# deciding which m to go next
				next_m = 0
				mm = m +1
				while mm < M and data[mm][0] == data[m][0]:
					# if also data[mm][9] == 0, it means a 2nd integron has been found, so break
					if data[mm][1] == "attC_site" and data[mm][9] == 0: 
						break
					if data[mm][1] == "attC_site":
						next_m = mm + 1
						attCind.append(mm)
					mm = mm + 1
				# end while
				# if no other attC 	
				n_remove = 0
				if not attCind:
					# if strand == '+' --> remove all
					if strand == '+' and m+1 <= mm:
						del data[(m+1):mm]
						remove_not_int += (mm - (m+1))
						mm = mm - (mm-(m+1))
						n_remove = mm-(m+1)
					# if strand == '-' --> remove all - 1 (-1 is orf in bottom cassette)
					elif strand == '-' and m+2 <= mm:
						del data[(m+2):mm]
						remove_not_int += (mm - (m+2))
						n_remove = mm-(m+2)
						mm = mm - (mm-(m+2))
					last = mm
				else:
					# correcting mm to be the 1st entry after the last attC on the integron.
					last = mm					
					mm = attCind[-1] + 1
					if strand == '+':
						# if CDSs exist after the last attC_site
						# --> remove all CDS after the last attCsite
						if attCind[-1] < (last-1):
							del data[(attCind[-1]+1):last]
							remove_not_int += (last - (attCind[-1]+1))
							mm = mm - (last-(attCind[-1]+1))
							last -= (last-(attCind[-1]+1))
					if strand == '-':
						# if 2+ CDSs exist after the last attC_site
						# --> remove all 2+ CDS after the last attCsite
						# here: mm is the next entry after the last attC, 
							#i.e. the CDS we want to keep!
						if attCind[-1] < (last-2):
							del data[(attCind[-1]+2):(last)]
							remove_not_int += (last - (attCind[-1]+2))
							mm = mm - (last-(attCind[-1]+2))
							last -= (last-(attCind[-1]+2))
				# next_m == 0, when only one attC
				if next_m == 0:
					m = mm
				else:
					m = last #next_m - n_remove
				M = len(data)
#		# ----- Remove attC singletons without CDS or CDS futher than d_CDS_attC ----- #
		# ----- Remove attC that do not have a CDS associated (empty cassettes) ------ #
		M = len(data)
		m = 0
		while m < M:
			# CASE 1: attC, PLUS
			# case 1A) m = 0, attC, +, and no previous CDS:
			if data[m][1] == "attC_site" and data[m][4] == "+" and m-1 < 0:
				del data[m]
				m -= 1
			# case 1B) attC_1, +, no previous CDS
			elif data[m][1] == "attC_site" and data[m][4] == "+" and m-1 >= 0 and data[m-1][0] != data[m][0]:
				del data[m]
				m -= 1
			# case 1C) attC, +, no previous CDS
			elif data[m][1] == "attC_site" and data[m][4] == "+" and m-1 >= 0 and data[m-1][1] != "CDS":
				del data[m]
				m -= 1
			# CASE 1: attC, MINUS
			# case 1A) m = M-1, attC, -:
			elif data[m][1] == "attC_site" and data[m][4] == "-" and m == M-1:
				del data[m]
				m -= 1
			# case 1B) attC_1, +, no next CDS
			elif data[m][1] == "attC_site" and data[m][4] == "-" and m+1 < M and data[m+1][0] != data[m][0]:
				del data[m]
				m -= 1
			# case 1C) attC, +, no next CDS
			elif data[m][1] == "attC_site" and data[m][4] == "-" and m+1 < M and data[m+1][1] != "CDS":
				del data[m]
				m -= 1
			M = len(data)
			m += 1 
		# -------------- adding info about CDS -------- #
		# each "integron" start in only 3 possible ways:
		# 1) CDS, attC, ...		(if plus strand)
		# 1b) CDS, CDS, attC	(if plus strand, and attC_1 was empty)
		# 2) attC, CDS, ... 	(if minus strand)
		# 2b) attC, CDS, CDS 	(if minus strand, and attC_1 was empty)
		# NOTE THAT: some attC_sites (in empty cassettes) have been removed, so, 
			#some single cassettes might look as singletons, but in reality there are not!
		M = len(data)
		m = 0
		while m < M:
			# case 1) attC, plus strand
			# 1+ CDS per attC_site, look for the previous attC_site, name all CDS_#
			#	where # is the attC_#
			if data[m][1] == "attC_site" and data[m][4] == "+":
				mm = m - 1
				while mm >= 0 and data[mm][1] != "attC_site":
					if data[mm][0] == data[m][0]:
						mm -= 1
					else:
						break
				for i in range(m-1,mm,-1):
					data[i][8] = "CDS_"+data[m][8].split("_")[1]
				m += 1
			# case 2) CDS, minus strand
			# 1+ CDS per attC_site, m-1 is always attC.
			# Look for the next CDSs, name all CDS_# where # is the attC_#
			elif data[m][1] == "attC_site" and data[m][4] == "-":
				mm = m + 1
				while mm < M and data[mm][1] != "attC_site":
					if data[mm][0] == data[m][0]:
						mm += 1
					else:
						break
				for i in range(m+1, mm):
					data[i][8] = "CDS_"+data[m][8].split("_")[1]
				m = mm
			else:
				m += 1
		# -------------- Completing with NA's --------- #
		M = len(data)
		for m in range(0,M):
			if data[m][1] == 'CDS':
				for i in range(12,19):
					data[m].insert(i, "NA")
		# -------------- recalculating distances ------ #
		names_CDS_complete = ""
		if len(data) > 0:
			data[0][7] = 'ini0'
			pos_attC = int(data[0][3])
			for m in range(1,M):
				if data[m-1][0] !=  data[m][0]:
					data[m][7] = 'ini0'
					if data[m][1] == "attC_site":
						pos_attC = int(data[m][3])
				elif data[m][1] == "attC_site" and data[m][2] - int(pos_attC) < dist_threshold:
					data[m][7] = data[m][2] - int(data[m-1][3])
					pos_attC = int(data[m][3])
				elif data[m][1] == "CDS" and data[m][2] - int(data[m-1][3]) < d_CDS_attC:
					data[m][7] = data[m][2] - int(data[m-1][3])
				else:
					data[m][7] = 'ini0'
		# -------------- saving ----------------------- #
		# name and ini_pos of ORFs extract info to save fasta later
		names = []
		names_attC = []
		fileout = open(out+".results",'w')
		fileout.write("id\telement\tstart\tend\tstrand\tscore\tlen\tdist\tarray_no\tdist_attC\tVscore\tR\'\tsp\'\tL\'\tloop\tL\"\tsp\"\tR\"\n")
		for m in range(0,M):
			if 'CDS' in data[m][8]:
				names.append(str(data[m][0])+"#"+str(data[m][2])+"#"+str(data[m][3])+"#"+str(data[m][4])+"#"+str(data[m][8]))
			elif 'attC' in data[m][8]:
				names_attC.append(str(data[m][0])+"#"+str(data[m][2])+"#"+str(data[m][3])+"#"+str(data[m][4])+"#"+str(data[m][8]))
			data[m][2] = str(data[m][2])
			data[m][6] = str(data[m][6])
			data[m][7] = str(data[m][7])
			data[m][9] = str(data[m][9])
			fileout.write("\t".join(data[m]))
			fileout.write("\n")
		fileout.close()
		# --------- Saving significant ORFs (FASTA) ----------- #
		if len(names) > 0:
			M = len(data)
			all_orf = SeqIO.parse(open("orfresults/"+out+"_ALLorf.faa"),'fasta')
			orf_fasta = "orfresults/"+out+"_orf.faa"
			with open(orf_fasta, 'w') as f:
				for seq in all_orf:
					sta = seq.description.split()[2]
					sto = seq.description.split()[4]
					acc = seq.description.split()[0]
					tap = seq.description.split()[6]
					if tap == '1':
						tap = '+'
					else:
						tap = '-'
					# PRODIGAL adds _# to each acc
					# removing it to match original				
					acc = acc.split("_")[:-1]
					acc = "_".join(acc)
					acc2 = acc+"#"+sta+"#"+sto+"#"+tap
					for n in range(0,len(names)):
						if acc2 in names[n]:
							seq.id = names[n].split()[0]
							seq.description = names[n]
							foo = SeqIO.write([seq],f,'fasta')
			f.close()
			#
			all_orf_nt = SeqIO.parse(open("orfresults/"+out+"_ALLorf.fna"),'fasta')
			orf_fasta_nt = "orfresults/"+out+"_orf.fna"
			with open(orf_fasta_nt, 'w') as fnt:
				for seq in all_orf_nt:
					sta = seq.description.split()[2]
					sto = seq.description.split()[4]
					acc = seq.description.split()[0]
					tap = seq.description.split()[6]
					if tap == '1':
						tap = '+'
					else:
						tap = '-'
					# PRODIGAL adds _# to each acc
					# removing it to match original				
					acc = acc.split("_")[:-1]
					acc = "_".join(acc)
					acc2 = acc+"#"+sta+"#"+sto+"#"+tap
					for n in range(0,len(names)):
						if acc2 in names[n]:
							seq.id = names[n].split()[0]
							seq.description = names[n]
							foo = SeqIO.write([seq],fnt,'fasta')
			fnt.close()
			#
			all_attC = SeqIO.parse(open("cmresults/"+out+"_infernal.fasta"),'fasta')
			attC_fasta = "cmresults/"+out+"_attC.fasta"
			with open(attC_fasta, 'w') as fattC:
				for seq in all_attC:
					acc = seq.description.split()[0]
					for n in range(0,len(names_attC)):
						tmp = names_attC[n].split("#")[0]
						if acc == tmp:
							name = names_attC[n]
							sta = int(name.split("#")[1])-1
							sto = int(name.split("#")[2])
							seq.description = name
							seq.id = name
							aux = seq.seq
							seq.seq = seq.seq[sta:sto]
							foo = SeqIO.write([seq],fattC,'fasta')
							seq.seq = aux
			fattC.close()
			#
			all_coord = open("cmresults/"+out+"_attC.res", 'rb')
			attC_coord = "cmresults/"+out+"_attC.coord"
			with open(attC_coord, 'w') as fcoord:
				for line in all_coord:
					tmp = line.split()
					acc = tmp[0]
					sta = tmp[1]
					sto = tmp[2]
					tap = tmp[3]
					acc = acc+'#'+sta+'#'+sto+'#'+tap+'#'
					for n in range(0,len(names_attC)):
						if acc in names_attC[n]:
							fcoord.write(line)
							break
			fcoord.close()

			os.unlink(output_file_faa)
			os.unlink(output_file_fna)

	# ------------------------------------------------------ #
	# ---------------------------------------------------- #
	# Saving info about removed ORFs

		with open(f"{output_dir}/{output_file}.filtering", 'a') as fileout:
			writer = csv.writer(fileout, delimiter="\t")
			writer.writerow(["Del_attC_ol:", str(remove_overlap_attC)])
			writer.writerow(["Total_orfs:", str(total_orf)])
			writer.writerow(["Del_orf_k:", str(remove_k)])
			writer.writerow(["Del_orf_ol:", str(remove_overlap)])
			writer.writerow(["Del_orf_ni:", str(remove_not_int)])
			writer.writerow(["ol= overlap, k= score_below_threshold, ni= not_integron, ws= wrong_strand"])	
	

#---------------------------------------------------------------------------#



#---------------------------------------------------------------------------#
def main(fastafile, output_directory, cm_model=None, both=True, nseq=1000, nthread=6, k_cm=20, k_orf=0, save_orf=True, dist_threshold=4000, d_CDS_attC=500):
	# creating outfile name	
	file_name = os.path.splitext(fastafile)[1]
	if file_name not in [".fasta", ".fa", ".fna"]:
		logging.error("Unknown file format, please input a .fasta, .fa or .fna file.")
		sys.exit(1)
	
	if not cm_model:
		my_resources = files('migfinder') / 'cm_model'
		cm_model = (my_resources / 'selection109_oriR.cm')

	# calling hattci
	hattci_fasta = hattci(fastafile, output_directory, nseq= nseq, nthread = nthread)
	
	logging.info("HattCI done!")

	# calling infernal + pos_processing
	infernal(hattci_fasta, output_directory, cm_model)

	logging.info("Infernal done!")

	posproc(fastafile, output_directory, k_cm)

	logging.info("Post proc done!")

	if os.path.exists("cmresults/" + out + "_attC.res"):
		prodigal(fastafile="../cmresults/" + out + "_infernal.fasta", out=out, save_orf=save_orf )
		logging.info("Prodigal done! Starting pos-processing...")
		posproc2(output_directory, k_orf)
		logging.info("Pos-processing done!")
	else:
		logging.info("No attC found, skipping Prodigal step...")
		

def migfinder_cli():
	parser = argparse.ArgumentParser(description="Metagenomic Integron-associated Gene finder")
	parser.add_argument("-f", "--fasta", required=True, help="Input fasta file. Valid extensions are .fasta, .fna, and .fa")
	parser.add_argument("-o", "--output", required=True, help="Output directory")	
	parser.add_argument("-c", "--cmmodel", required=False, help="Covariance model used by Infernal to validate the attC site secondary structure [default=None]")
	parser.add_argument("-b", "--strand", required=False, help="Perform HattCI in both strands [default=True]")
	parser.add_argument("-n", "--seqs", required=False, help="Number of sequences processed at a time by HattCI [default=1000]")
	parser.add_argument("-t", "--threads", required=False, help="Number of threads to run HattCI, [default=6]")	
	parser.add_argument("-e", "--score", required=False, help="Threshold used to filter Infernal results [default=20]")
	parser.add_argument("-r", "--orf", required=False, help="Threshold used to filter HattCI results [default=0]")
	parser.add_argument("-s", "--save", required=False, help="Save ORF results in a separate fasta file [default=True]")
	parser.add_argument("-a", "--adist", required=False, help="Max distance allowed to consider two adjacent attC sites part of the same integron [default=4000]")	
	parser.add_argument("-d", "--odist", required=False, help="Max distance allowed between ORF and attC site in the same gene cassette [default=500]")
	
	args = parser.parse_args()

	return main(args.fasta, args.output, args.cmmodel, args.strand , args.seqs , args.threads , args.score , args.orf , args.save , args.adist , args.odist)


if __name__ == "__main__":
	migfinder_cli()
#---------------------------------------------------------------------------#
