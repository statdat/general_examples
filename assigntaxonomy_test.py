
'''
testing framework for Taxonomic Assignment
'''

import os
from os.path import join as pjoin

from timeit import default_timer as timer

import unittest

os.chdir('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/scripts/enr_analysis')



base_test_dir ='/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/assign_taxonomy'



#####-----------------------------------test 1 setup #####------------------------------------------------------
#creating inputs for extract_region_attributes.py
input_data = pjoin(base_test_dir, 'test_files')
base_dir = pjoin(base_test_dir, 'test1')
os.system('rm -R %s' % base_dir)
os.system('mkdir %s' % base_dir)

#####-----------------------------------test 1 ExtractMarkers.py #####------------------------------------------------------
from ExtractMarkers import *
input_marker_query_file = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/assign_taxonomy/test_files/gyrb_rpob_rpod.fasta'
processes = 8
run_ExtractMarkers(INPUT_DATA_PATH_=input_data,DATA_PATH_=base_dir,INPUT_MARKER_QUERY_FILE_=input_marker_query_file,processes=processes)


#####-----------------------------------test 2 set up #####------------------------------------------------------
test1_workdir = base_dir
base_dir = pjoin(base_test_dir, 'test2')
os.system('rm -R %s' % base_dir)
os.system('mkdir %s' % base_dir)
os.system('cp -a %s/. %s' % (test1_workdir, base_dir))

#####-----------------------------------test 2 MatchSpecies.py #####------------------------------------------------------
from MatchSpecies import *
processes = 8
type_strains_df = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/assign_taxonomy/test_files/type_strains.csv'
run_MatchSpecies(DATA_PATH_=base_dir,INPUT_TYPE_SPECIES_DF_=type_strains_df,processes=processes)



#####-----------------------------------test 3 set up #####------------------------------------------------------
test1_workdir = base_dir
base_dir = pjoin(base_test_dir, 'test3')
os.system('rm -R %s' % base_dir)
os.system('mkdir %s' % base_dir)
os.system('cp -a %s/. %s' % (test1_workdir, base_dir))

#####-----------------------------------test 3 MSAMarkers.py #####------------------------------------------------------
from MSAMarkers import *
processes = 8
input_cf_df = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/assign_taxonomy/test_files/processed_genomes_metadata.csv'
input_delim = 't1'
run_MSAMarkers(DATA_PATH_=base_dir, INPUT_CF_DF_=input_cf_df,INPUT_DELIMITER=input_delim,processes=processes)










