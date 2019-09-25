# This script configures the directory and files needed for running the mhc step.

import os
from optparse import OptionParser
from shutil import copyfile

# Parsing arguments from command line
parser = OptionParser(usage='python configuration.py <options>')
parser.add_option('--directory', dest='directory', action='store')
parser.add_option('--patientID', dest='patientID', action='store')
parser.add_option('--peptides_directory', dest='peptides_directory', action='store')
(options, args) = parser.parse_args()

# Assign arguments
directory = options.directory
peptides_dir = os.path.join(directory, 'peptides')
hla_dir = os.path.join(directory, 'hla')
iebd_dir = os.path.join(directory, 'IEDB_out')

# Make directory
if not os.path.exists(directory):
    os.mkdir(directory)

if not os.path.exists(peptides_dir):
    os.mkdir(peptides_dir)

if not os.path.exists(hla_dir):
    os.mkdir(hla_dir)

if not os.path.exists(iebd_dir):
    os.mkdir(iebd_dir)

# Copy peptide files
peptide_lengths = [15, 17, 19, 21]
for i in peptide_lengths:
    original_peptide_file = os.path.join(options.peptides_directory, options.patientID + '.VarScan_vep.' + str(i) + '.peptide')
    new_peptide_file = os.path.join(peptides_dir, options.patientID + '.' + str(i) + '.txt')
    copyfile(src=original_peptide_file, dst=new_peptide_file)