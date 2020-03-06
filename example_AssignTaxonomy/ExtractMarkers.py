
'''
INPUT_DATA_PATH is the base directory containing the prokka and blast files produced by fabv_InputCreation.py
DATA_PATH is the output directory
USER must stupply a file containing a fasta seq of gyrb, rpob, and rpod.
Outputs: alll gyrb,rpob,rpod sequences concatenated indivudually and a long contaentated gyrb,rpob,rpod sequence
per cf input. 
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



def get_AssignTaxonomyPaths(INPUT_DATA_PATH_,DATA_PATH_,INPUT_MARKER_QUERY_FILE_):
	'''
	'''
# INPUT_DATA_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test_files/'
# DATA_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test1'
# INPUT_MARKER_QUERY_FILE = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test_files/gyrb_rpob_rpod.fasta'
	global INPUT_DATA_PATH, DATA_PATH, INPUT_MARKER_QUERY_FILE
	INPUT_DATA_PATH = INPUT_DATA_PATH_
	DATA_PATH = DATA_PATH_
	INPUT_MARKER_QUERY_FILE = INPUT_MARKER_QUERY_FILE_


def parallel_ExtractMarkerSeqs(cf,INPUT_DATA_PATH,DATA_PATH,INPUT_MARKER_QUERY_FILE,hits_output,concatseq_output):
	'''
	##Step1: blast
	##Step2: make a hit table of blast results and filter for 95% coverage and three unique hits. Get the highest bit score from each
	##Step3: Extract nucleotide sequences of the best hits from the hit table
	##Step4: Write the query-specific nucleotide sequences to corresponding output directory
	##Step5: Write the concatenated sequence for the cf to the corresponding concatseq_output directory
	'''
	##Step1: blast for inputted query seqs
	cf_blast_db = pjoin(INPUT_DATA_PATH,'blast_output',cf) #location of blast database files
	cf_hits_output = pjoin(DATA_PATH,hits_output,cf)
	os.system('blastp -query %s -db %s -max_target_seqs 100 -evalue 1e-6 -outfmt "10 sseqid qseqid mismatch positive gaps ppos pident qcovs evalue bitscore sstart send" -num_threads 1  -out %s.csv' % (INPUT_MARKER_QUERY_FILE, cf_blast_db, cf_hits_output))
	##Step2: Read in hit table and check for hits meeting criteria: 1 of each, and 95% coverage. Save the reduced df
	col_names = ['sseqid', 'qseqid', 'mismatch', 'positive', 'gaps', 'ppos', 'pident', 'qcovs', 'evalue', 'bitscore', 'sstart', 'send']
	df_cf = pd.read_csv(cf_hits_output+'.csv',names=col_names)
	#filter cf_df for >=95% coverage
	df_cf = df_cf[df_cf['qcovs']>=95]
	#check for 3 hits. if there are not three hits then delete the file and exit the function
	if len(df_cf['qseqid'].unique())!=3:
		os.system('rm %s.csv'%cf_hits_output)
		print(cf+' failed to find all three genes')
		return None
	#extract the max bitsore from each unique group
	#the lambda x is basically an 'i' in a for loop but for a dataframe based on the groupby
	df_cf = df_cf.groupby('qseqid').apply(lambda x: x[x['bitscore'] == x['bitscore'].max()])
	#dropping duplicates: this is for cases when two sseqids from the same qseqid have equal bitscores
	df_cf = df_cf.drop_duplicates(['qseqid'])
	#adding cleaned_filename for future reference
	df_cf['cleaned_filename'] = cf
	#numbering dataframe as gene_order to specify concatenation later
	df_cf['gene_order'] = list(range(1,4))
	#save df_cf with same name as before
	df_cf.to_csv(cf_hits_output+'.csv',index=False,header=False)
	##Step3: extracting nucleotide sequence to a database containing sseqid, qseqid, and gene order
	df_query_seq = []
	for record in SeqIO.parse(pjoin(INPUT_DATA_PATH,'prokka_output',cf+'.ffn'),'fasta'):
		if record.id in df_cf['sseqid'].tolist():
			#extract the query name
			query_id = df_cf[df_cf['sseqid']==record.id]['qseqid'].item()
			sseqid_seq = str(record.seq)
			#the sseqid(record.id), query_id, the sseqid sequence
			df_query_seq.append([record.id, query_id, sseqid_seq])
	#converting df_query_seq to array and then pandas df
	df_query_seq = np.asarray(df_query_seq)
	df_query_seq = pd.DataFrame(df_query_seq)
	df_query_seq.columns = ['sseqid','qseqid','sequence']
	##merging the gene order from df_cf with df_query_seq
	df_cf_order = df_cf[['qseqid','gene_order']].reset_index(drop = True) #dropping is necessary
	#df_query_seq is sorted by gene_order 1 - 3
	df_query_seq = df_query_seq.merge(df_cf_order, on = 'qseqid').sort_values(by=['gene_order'],ascending=True)
	##Step4: Writing to each sequence to its respective folder location
	for sseqid in df_query_seq['sseqid'].tolist():
		#sequence
		sseqid_seq = df_query_seq[df_query_seq['sseqid']==sseqid]['sequence'].item()
		#query id
		query_id = df_query_seq[df_query_seq['sseqid']==sseqid]['qseqid'].item()
		#file output name
		query_output = pjoin(DATA_PATH,query_id+'_output',sseqid)
		with open(query_output+'.txt','w') as query_outfile:
			#the fasta file has the sseqid as record.id, cleaned filename (cf) as record.descrption, and record.seq is sequence
			query_outfile.write('>%s %s\n%s\n'%(sseqid,cf,sseqid_seq))
	##Step5: Writing a concatenated sequence of all sseqid in the order of gene_order
	concatenated_seq = ''
	for seq in df_query_seq['sequence'].tolist():
		concatenated_seq += seq
	#concatenated seq being written to concatseq_output
	with open(pjoin(concatseq_output,cf)+'.txt','w') as concatseq_outfile:
		#record.id is cleaned_filename (cf). record.seq is concatenated sequence
		concatseq_outfile.write('>%s\n%s\n'%(cf,concatenated_seq))



def extract_MarkerSeqs(DATA_PATH,INPUT_DATA_PATH,INPUT_MARKER_QUERY_FILE,processes):
	'''
	Step1:make output directories Step2: make parallel commands
	'''
	##Making output directories
	hits_output = pjoin(DATA_PATH,'hits_output')
	concatseq_output = pjoin(DATA_PATH,'concatseq_output')
	#these output folders have the same names as the defline in INPUT_MARKER_QUERY_FILE+'_output'
	gyrb_output = pjoin(DATA_PATH,'pao1_gyrb_output')
	rpob_output = pjoin(DATA_PATH,'pao1_rpob_output')
	rpod_output = pjoin(DATA_PATH,'pao1_rpod_output')
	output_directories = [hits_output,concatseq_output,gyrb_output,rpob_output,rpod_output]
	[os.mkdir(x) for x in output_directories if os.path.exists(x)==False]
	##Making parallel commands
	#all cf will be submitted
	df_allcf = pd.read_csv(pjoin(INPUT_DATA_PATH,'processed_genomes_metadata.csv'))
	#parallel commands
	parallel_input=[[cf,INPUT_DATA_PATH,DATA_PATH,INPUT_MARKER_QUERY_FILE,hits_output,concatseq_output] for cf in df_allcf['cleaned_filename'].tolist()[:]]
	pool = Pool(processes = processes)
	pool.starmap(parallel_ExtractMarkerSeqs,parallel_input)
	pool.close()


def concatenate_ExtractMarkerSeqsOutputs(DATA_PATH):
	'''
	Step1:All hits files will be concatenated and saved into the main DATA_PATH folder and the rest of the files + directory deleted
	Step2:All files in query-specific directories can be concatenated into one large file and saved into the main DATA_PATH folder and the rest of the files deleted + directory deleted
	'''
	##step1:
	hits_output = pjoin(DATA_PATH,'hits_output')
	#hit files to concatenate
	concat_hits_files = [pjoin(hits_output,x) for x in os.listdir(hits_output) if x.endswith('.csv')]
	#hit output filename
	concat_hits_output = pjoin(DATA_PATH,'concatenated_hits')
	cat_Files(concat_hits_files,'csv',DATA_PATH,concat_hits_output)
	#saving hit output
	concat_hits_cols = ['sseqid', 'qseqid', 'mismatch', 'positive', 'gaps', 'ppos', 'pident', 'qcovs', 'evalue', 'bitscore', 'sstart', 'send','cleaned_filename','gene_order']
	df_concat_hits = pd.read_csv(concat_hits_output+'.csv',names=concat_hits_cols)
	df_concat_hits.to_csv(concat_hits_output+'.csv',index=False)
	#delete hts directory
	os.system('rm -R %s'%(hits_output))
	##Step2
	#directories to concatenate
	gyrb_output = pjoin(DATA_PATH,'pao1_gyrb_output')
	rpob_output = pjoin(DATA_PATH,'pao1_rpob_output')
	rpod_output = pjoin(DATA_PATH,'pao1_rpod_output')
	#looping through qseqid ids
	for qseqid in ['pao1_gyrb','pao1_rpob','pao1_rpod']:
		#query files to concatenate
		all_qseqid_files = [pjoin(DATA_PATH,qseqid+'_output',x) for x in os.listdir(pjoin(DATA_PATH,qseqid+'_output')) if x.endswith('.txt')]
		#query output filename
		qseqid_output_file = pjoin(DATA_PATH,qseqid)
		cat_Files(all_qseqid_files,'txt',DATA_PATH,qseqid_output_file)
	#deleting files
	os.system('rm -R %s %s %s'%(gyrb_output,rpob_output,rpod_output))



def run_ExtractMarkers(INPUT_DATA_PATH_,DATA_PATH_,INPUT_MARKER_QUERY_FILE_,processes):
	s = timer()
	get_AssignTaxonomyPaths(INPUT_DATA_PATH_,DATA_PATH_,INPUT_MARKER_QUERY_FILE_)
	print('extract_MarkerSeqs')
	extract_MarkerSeqs(DATA_PATH,INPUT_DATA_PATH,INPUT_MARKER_QUERY_FILE,processes)
	print('concatenate_ExtractMarkerSeqsOutputs')
	concatenate_ExtractMarkerSeqsOutputs(DATA_PATH)
	e = timer()
	print('total time:')
	print(e-s)


















