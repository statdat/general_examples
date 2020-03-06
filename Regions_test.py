
'''
testing framework for isolate_fab_homolgs.py
'''

import os
from os.path import join as pjoin

from timeit import default_timer as timer

import unittest

os.chdir('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/scripts/enr_analysis')



base_test_dir ='/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_regs'



#####-----------------------------------test 1 setup #####------------------------------------------------------



#creating inputs for extract_region_attributes.py
input_data = pjoin(base_test_dir, 'test_files')
base_dir = pjoin(base_test_dir, 'test1')
query_df_dir = pjoin(base_test_dir, 'test_files')
os.system('rm -R %s' % base_dir)
os.system('mkdir %s' % base_dir)

# #####-----------------------------------test 1 fabv_InputCreation.py #####------------------------------------------------------
from fabv_InputCreation import *
processes = 8
run_InputCreation(INPUT_DATA_PATH_=input_data,DATA_PATH_=base_dir,INPUT_FILE_DF_=pjoin(query_df_dir,'isolate_info.csv'),processes_=processes)



# #####-----------------------------------test 2 set up #####------------------------------------------------------
# test1_workdir = base_dir
# base_dir = pjoin(base_test_dir, 'test2')
# out_dir = base_dir
# input_query_fasta_file = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/tests/fabv_regs/test_files/pao1_fabv.fasta'
# os.system('rm -R %s' % base_dir)
# os.system('mkdir %s' % base_dir)
# os.system('cp -a %s/. %s' % (test1_workdir, base_dir))

# #####-----------------------------------test 2 fabv_GetHomologs.py #####------------------------------------------------------
# from fabv_GetHomologs import *
# processes=8
# input_query_name = 'pao1_fabv'

# run_GetHomologs(INPUT_DATA_PATH_=base_dir,DATA_PATH_=base_dir,INPUT_QUERY_FASTA_FILE_=input_query_fasta_file,INPUT_QUERY_NAME_=input_query_name,processes_=processes)

