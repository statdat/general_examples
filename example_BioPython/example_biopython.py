#!//anaconda/bin/python

"""
Make sure the above has been changed for your own path

"""


import os
from Bio import SeqIO


for record in SeqIO.parse("example.fasta", "fasta"):
	#this prints the header id of the fasta
	print (record.id)


for record in SeqIO.parse("example.fasta", "fasta"):
	#this prints the description of the fasta 
	print (record.description)


for record in SeqIO.parse("example.fasta", "fasta"):
	#this prints the sequence of the fasta
	print (record.seq)



#opening a fasta file and writing everything back originally
with open("processed_example_1.fasta", "w") as outfile:
	for record in SeqIO.parse("example.fasta", "fasta"):
		description = record.description
		sequence = record.seq
		#description and sequence need to explicitly be converted to a str object using str()
		outfile.write('>'+str(description)+'\n'+str(sequence)+'\n')
outfile.close()

#opening a fasta file and writing a new description
with open("processed_example_2.fasta", "w") as outfile:
	for record in SeqIO.parse("example.fasta", "fasta"):
		description = "my_new_description"
		sequence = record.seq
		#description and sequence need to explicitly be converted to a str object using str()
		outfile.write('>'+str(description)+'\n'+str(sequence)+'\n')
outfile.close()


#opening a fasta file and writing a new description and only the first 10 nucleotides
with open("processed_example_3.fasta", "w") as outfile:
	for record in SeqIO.parse("example.fasta", "fasta"):
		description = "my_new_description"
		sequence = record.seq
		sequence_string = str(sequence)
		sequence_string_short = sequence_string[:10]
		#description and sequence need to explicitly be converted to a str object using str()
		outfile.write('>'+str(description)+'\n'+str(sequence_string_short)+'\n')
outfile.close()













