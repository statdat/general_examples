
'''
DATA_PATH is the input and output directory. uses the files outputted by fabv_ExtractMarkers.py
USER must supply a dataframe with the type species used for taxonomic reassignment process.
Output is a table containing reassigned taxonomies: taxonomic_assignments.csv
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



def get_MatchSpeciesPATHS(DATA_PATH_,INPUT_TYPE_SPECIES_DF_):
	'''
	'''
# DATA_PATH = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test2'
# INPUT_TYPE_SPECIES_DF = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_asstax/test_files/type_strains.csv'
	global DATA_PATH, INPUT_TYPE_SPECIES_DF
	DATA_PATH = DATA_PATH_
	INPUT_TYPE_SPECIES_DF = INPUT_TYPE_SPECIES_DF_


def make_BlastReferenceDb(DATA_PATH, INPUT_TYPE_SPECIES_DF):
	'''
	Create a reference database of concatenated sequences belonging to type species
	'''
	blast_reference_dir = pjoin(DATA_PATH,'blast_reference')
	if os.path.exists(blast_reference_dir)==False: os.mkdir(blast_reference_dir)
	concatseq_input_dir = pjoin(DATA_PATH,'concatseq_output')
	##concatenate all type species concatenate sequences into one large fasta file called type_species_reference.txt
	df_type_species = pd.read_csv(INPUT_TYPE_SPECIES_DF)
	#adding .txt to the end of all filenames so that it can be searched in concatseq_input_dir in the subsequent line
	type_species_list = [x+'.txt' for x in df_type_species['cleaned_filename']]
	#files to concatenate:
	type_species_file_list = [pjoin(concatseq_input_dir,x) for x in os.listdir(concatseq_input_dir) if x in type_species_list]
	#output file name:
	type_species_output = pjoin(DATA_PATH,'type_species_reference')
	cat_Files(type_species_file_list,'txt',DATA_PATH,type_species_output)
	##make blast reference database
	#name for blast reference database
	type_species_blast = pjoin(blast_reference_dir,'type_species_reference')
	os.system('makeblastdb -in %s.txt -out %s -dbtype nucl -title "type_species_reference_db" -parse_seqids' %(type_species_output,type_species_blast))



def parallel_MatchBestSpeciesHit(cf,DATA_PATH,concatseq_hits_dir,type_species_blast_db):
	'''
	Step1:Blast, step2:filter hits table for >=95 qcovs and max bitscore and then save file
	'''
	##Step1:blast
	cf_blast_input = pjoin(DATA_PATH,'concatseq_output',cf)
	cf_blast_output = pjoin(concatseq_hits_dir,cf)
	os.system('blastn -query %s.txt -db %s -evalue 1e-8 -outfmt "10 sseqid qseqid mismatch positive gaps ppos pident qcovs evalue bitscore sstart send" -num_threads 1  -out %s.csv' % (cf_blast_input, type_species_blast_db, cf_blast_output))
	##Step2:reduce dataframe to best hit that is also 95% or greater qcovs
	col_names = ['sseqid', 'qseqid', 'mismatch', 'positive', 'gaps', 'ppos', 'pident', 'qcovs', 'evalue', 'bitscore', 'sstart', 'send']
	df_cf = pd.read_csv(cf_blast_output+'.csv',names=col_names)
	#filter cf_df for >=95% coverage
	df_cf = df_cf[df_cf['qcovs']>=95]
	#filter for greatest bitscore
	df_cf = df_cf[df_cf['bitscore']==df_cf['bitscore'].max()]
	#dropping duplicates: this is for cases when two sseqids from the same qseqid have equal bitscores
	df_cf = df_cf.drop_duplicates(['qseqid'])
	#saving
	df_cf.to_csv(cf_blast_output+'.csv',index=False,header=False)



def match_BestSpeciesHit(DATA_PATH,processes):
	'''
	Blast each individual concatenated sequence against type_species_reference blast database
	'''
	#making output directory with the concat sequence blast hits
	concatseq_hits_dir = pjoin(DATA_PATH,'concatseq_hits')
	if os.path.exists(concatseq_hits_dir)==False: os.mkdir(concatseq_hits_dir)
	#concatenated_hits.csv contains blast info for all cf that passed the qcov and qseqid number thresholds
	df_cf = pd.read_csv(pjoin(DATA_PATH,'concatenated_hits.csv'))
	#reducing to unique cfs
	cf_list = df_cf['cleaned_filename'].unique().tolist()
	#the reference database filename
	type_species_blast_db = pjoin(DATA_PATH,'blast_reference','type_species_reference')
	##Parallel commands
	parallel_input = [[cf,DATA_PATH,concatseq_hits_dir,type_species_blast_db] for cf in cf_list]
	pool = Pool(processes = processes)
	pool.starmap(parallel_MatchBestSpeciesHit, parallel_input[:])
	pool.close()


def concatenate_MatchBestSpeciesHit(DATA_PATH):
	'''
	concatenate files in the concatseq_hits folder produced by match_BestSpeciesHits
	Step1: concatenate files Step2: delete .csv hits files folder and the .txt concat seqs files folder
	'''
	##Step1
	#list of all .csv output files
	all_concatseq_hit_files = [pjoin(DATA_PATH,'concatseq_hits',x) for x in os.listdir(pjoin(DATA_PATH,'concatseq_hits')) if x.endswith('.csv')]
	#output filename
	concatenated_concatseq_hits = pjoin(DATA_PATH,'taxonomic_assignments')
	cat_Files(all_concatseq_hit_files,'csv',DATA_PATH,concatenated_concatseq_hits)
	#adding headers to taxonomic_assignments.csv and renaming qseqid as cleaned_fileanme, sseqid_id as match
	col_names = ['match', 'cleaned_filename', 'mismatch', 'positive', 'gaps', 'ppos', 'pident', 'qcovs', 'evalue', 'bitscore', 'sstart', 'send']
	df_taxass = pd.read_csv(concatenated_concatseq_hits+'.csv',names=col_names)
	df_taxass.to_csv(concatenated_concatseq_hits+'.csv',index=False)
	##Step2: remove uneccesary folders: blast reference for searching, the concatenated sequences, the hits tables of concatenated sequnece blast
	r1 = pjoin(DATA_PATH,'blast_reference')
	r2 = pjoin(DATA_PATH,'concatseq_hits')
	r3 = pjoin(DATA_PATH,'concatseq_output')
	os.system('rm -R %s %s %s'%(r1,r2,r3))


def reassign_Taxonomy(DATA_PATH, INPUT_TYPE_SPECIES_DF):
	'''
	Takes the species name from INPUT_TYPE_SPECIES_DF and assigns it to taxonomic_assignments.csv
	'''
	##type species df: getting only cleaned_filename,type,and species
	df_input_type = pd.read_csv(INPUT_TYPE_SPECIES_DF)
	df_input_type = df_input_type[['cleaned_filename','species']].reset_index(drop = True)
	#renaming cleaned_filename to 'match' for later merging
	df_input_type = df_input_type.rename({'cleaned_filename':'match'},axis=1)
	#taxonomic assignments df
	df_taxass = pd.read_csv(pjoin(DATA_PATH,'taxonomic_assignments.csv'))
	#left merge on df_taxass
	df_taxass = df_taxass.merge(df_input_type,on='match')
	#saving dataframe
	df_taxass.to_csv(pjoin(DATA_PATH,'taxonomic_assignments.csv'),index=False)


def run_MatchSpecies(DATA_PATH_,INPUT_TYPE_SPECIES_DF_,processes):
	s = timer()
	get_MatchSpeciesPATHS(DATA_PATH_,INPUT_TYPE_SPECIES_DF_)
	print('make_BlastReferenceDb')
	make_BlastReferenceDb(DATA_PATH, INPUT_TYPE_SPECIES_DF)
	print('match_BestSpeciesHit')
	match_BestSpeciesHit(DATA_PATH,processes)
	print('concatenate_MatchBestSpeciesHit')
	concatenate_MatchBestSpeciesHit(DATA_PATH)
	print('reassign_Taxonomy')
	reassign_Taxonomy(DATA_PATH, INPUT_TYPE_SPECIES_DF)
	e = timer()
	print('total time:')
	print(e-s)









