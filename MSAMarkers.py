
'''
DATA_PATH is the input and output directory. Uses files produced by fabv_ExtractMarkers.py
USER must supply an output file delimiter
USER must supply a dataframe containing the cf that will be used to create a mafft alignment of the 
gyrb,rpob, and rpod sequences.
Output is an alignment of the above sequences that can then be used for tree building.
'''

import os
from os.path import join as pjoin

import re

from multiprocessing import Pool
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

import pandas as pd 
import numpy as np

import textwrap
from timeit import default_timer as timer


def cat_Files(file_list, extension, DATA_PATH_, out_filename_path):
	'''
	General purpose function to cat thousands of files and delete the intermediate chunking files
	file_list must have the path joined to the name and also the file extension
	extension should be written without the <.> examples: <txt> <csv>
	DATA_PATH_ is the output directory for the chunks
	out_filename_path must have the path joined to the name but not the extension
	'''
	major_length = len(file_list) - len(file_list)%100  #get length to chunk that is divisble by 100
	i = 0
	while i < len(file_list):
		i2 = i
		i += 100
		if i == major_length:
			print('final', '   ', len(file_list[i2:])) #for tracking
			chunk = ' '.join(file_list[i2:])
			chunk_path = pjoin(DATA_PATH_, 'final_chunk')
			os.system('cat %s > %s.%s' % (chunk, chunk_path, extension))
			break
		chunk = ' '.join(file_list[i2:i])
		print(i, '   ', len(file_list[i2:i])) #for tracking
		chunk_path = pjoin(DATA_PATH_, '%s_chunk' % str(i))
		os.system('cat %s > %s.%s' % (chunk, chunk_path, extension))
	wild_chunk_path = pjoin(DATA_PATH_,'*_chunk.%s' % extension)
	os.system('cat %s > %s.%s' % (wild_chunk_path, out_filename_path, extension))
	os.system('rm %s' % wild_chunk_path)


def get_MSAMarkersPATHS(DATA_PATH_, INPUT_CF_DF_):
	'''
	'''
# DATA_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test3'
# INPUT_CF_DF = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test_files/processed_genomes_metadata.csv'
	global DATA_PATH, INPUT_CF_DF
	DATA_PATH = DATA_PATH_
	INPUT_CF_DF = INPUT_CF_DF_


def make_MasterDf(DATA_PATH,INPUT_CF_DF):
	'''
	Filters the concatenated_hits df with cfs that passed all taxonomic criteria with the desired cf from the input list
	This df is then used for the rest of the script and is referred to as df_master
	'''
	##reading in df containing the desired cf inputs
	df_cf = pd.read_csv(INPUT_CF_DF)
	cf_input_list = df_cf['cleaned_filename'].tolist()[:]
	##reading in the master df contaiing cf that were able to have all three genes identified
	df_master = pd.read_csv(pjoin(DATA_PATH,'concatenated_hits.csv'))
	#filtering df_master for the desired cfs 
	df_master = df_master[df_master['cleaned_filename'].isin(cf_input_list)]
	#returning df_master that will be used by subsequent functions
	return(df_master)


def concatenate_SelectCf(df_master, DATA_PATH, INPUT_DELIMITER):
	'''
	Only files matching the selected cf will be concatenated. Supplied df only needs 'cleaned_filename' column to specify
	No selection for specific cf is done here or in this script
	'''
#INPUT_DELIMITER = 't1'
	##
	for qseqid in ['pao1_gyrb','pao1_rpob','pao1_rpod']:
		#reducing df_qseqid to df_specifc which only contains that matching qseqid 
		df_specific = df_master[df_master['qseqid']==qseqid]
		##Using df_specific to get all sseqids which match the record.ids in the specific merged txt files
		#writing to a new file which is unique due to the user INPUT_DELIMITER
		query_outfile = pjoin(DATA_PATH,INPUT_DELIMITER+'_'+qseqid)
		with open(query_outfile+'.txt','w') as outfile:
			#list of sseqids to search for
			sseqid_list = df_specific['sseqid'].tolist()
			for record in SeqIO.parse(pjoin(DATA_PATH,qseqid+'.txt'),'fasta'):
				if record.id in sseqid_list:
					outfile.write('>%s\n%s\n'%(record.id,str(record.seq)))


def align_IndividualSelectCf(INPUT_DELIMITER,processes):
	'''
	Each qseqid is now aligned 
	'''
	for qseqid in ['pao1_gyrb','pao1_rpob','pao1_rpod']:
		print('aligning %s'%qseqid)
		qseqid = INPUT_DELIMITER+'_'+qseqid
		qseqid_path = pjoin(DATA_PATH,qseqid)
		os.system('mafft --thread %s %s.txt > %s.aln'%(processes,qseqid_path,qseqid_path))


def parallel_ExtractCfAlignments(cf,temp_aln_output,DATA_PATH,INPUT_DELIMITER,df_master):
	'''
	'''
	##subset df_masfer to only rows with cf
	df_cf = df_master[df_master['cleaned_filename']==cf]
	#because df_master is already ordered in the correct contenation order, the subsetted df is already ordered
	#therefore, can loop through the sseqids and be sure they are in the correct concatenation order
	#adding aligned sequences to conc_aln_seq 
	conc_aln_seq = ''
	for sseqid in df_cf['sseqid'].tolist():
		#getting the name of the qseqid for which the aln file will be opened
		qseqid = df_cf[df_cf['sseqid']==sseqid]['qseqid'].item()
		qseqid_path = pjoin(DATA_PATH,INPUT_DELIMITER+'_'+qseqid)
		for record in SeqIO.parse(qseqid_path+'.aln','fasta'):
			if record.id == sseqid:
				conc_aln_seq += str(record.seq)
				break #exit loop for speed
	##Writing sequence to file
	#concat_aln_seq needs to have text be wrapped or else the hyphens will create line breaks
	conc_aln_seq = textwrap.wrap(conc_aln_seq,60,break_on_hyphens=False)
	with open(pjoin(temp_aln_output,cf+'.aln'),'w') as outfile:
		outfile.write('>%s\n'%cf)
		[outfile.write('%s\n'%i) for i in conc_aln_seq]


def extract_CfAlignments(df_master,DATA_PATH,INPUT_DELIMITER,processes):
	'''
	For each cf extract the alignments from each of the qseqid file outputs and contanenates the alignments into a single massive alignment
	'''
	##output folder for all the cf alignments
	temp_aln_output = pjoin(DATA_PATH,'temp_aln_output')
	if os.path.exists(temp_aln_output) == True: os.system('rm -R %s'%temp_aln_output)
	os.mkdir(temp_aln_output)
	##Making parallel commands
	parallel_input = [[cf,temp_aln_output,DATA_PATH,INPUT_DELIMITER,df_master] for cf in df_master['cleaned_filename'].tolist()]
	pool = Pool(processes = processes)
	pool.starmap(parallel_ExtractCfAlignments, parallel_input[:])
	pool.close()


def concatenate_CfAlignments(DATA_PATH, INPUT_DELIMITER):
	'''
	##Step1:concatenates the aln files produced by extract_CfAlignments and saves them in DATA_PATH
	##Ste2:removes the temp_aln_output directory
	'''
	##Step1
	#input file list
	input_aln_file_list = [pjoin(DATA_PATH,'temp_aln_output',x) for x in os.listdir(pjoin(DATA_PATH,'temp_aln_output')) if x.endswith('.aln')]
	#output filename
	conc_aln_output = pjoin(DATA_PATH,'mlsa_'+INPUT_DELIMITER)
	cat_Files(input_aln_file_list,'aln',DATA_PATH,conc_aln_output)
	##Step2
	os.system('rm -R %s'%(pjoin(DATA_PATH,'temp_aln_output')))



def run_MSAMarkers(DATA_PATH_, INPUT_CF_DF_,INPUT_DELIMITER,processes):
	s = timer()
	get_MSAMarkersPATHS(DATA_PATH_, INPUT_CF_DF_)
	print('make_MasterDf')
	df_master = make_MasterDf(DATA_PATH,INPUT_CF_DF)
	print('concatenate_SelectCf')
	concatenate_SelectCf(df_master, DATA_PATH, INPUT_DELIMITER)
	print('align_IndividualSelectCf')
	align_IndividualSelectCf(INPUT_DELIMITER,processes)
	print('extract_CfAlignments')
	extract_CfAlignments(df_master,DATA_PATH,INPUT_DELIMITER,processes)
	print('concatenate_CfAlignments')
	concatenate_CfAlignments(DATA_PATH, INPUT_DELIMITER)
	e = timer()
	print('total time:')
	print(e-s)















