# In this script, neoepitopes for each hla type are filtered:
# 1. Remove duplicates. These are not technically duplicates because the transcript names are different.
# However, because the gene names are the same, we will just keep one rows.
# 2. Only keep the neoepitopes where the score is less than 500

from optparse import OptionParser
import os

# Parsing arguments from command line
parser = OptionParser(usage='python filter_neoepitope.py <options>')
parser.add_option('--patientID', dest='patientID', action='store')
parser.add_option('--hla_type', dest='hla_type', action='store')
(options, args) = parser.parse_args()

input_fn = os.path.join(options.patientID, 'IEDB_out', options.hla_type, 'output_IEDB.15.17.19.21_minscore_minpeptide.txt')
rmdups_output_fn = os.path.join(options.patientID, 'IEDB_out', options.hla_type, 'output_IEDB.15.17.19.21_minscore_minpeptide_rmdups.txt')
scorefiltered_output_fn = os.path.join(options.patientID, 'IEDB_out', options.hla_type, 'output_IEDB.15.17.19.21_minscore_minpeptide_rmdups_scorefiltered.txt')

rmdups_outfile = open(rmdups_output_fn, 'w')
scorefiltered_outfile = open(scorefiltered_output_fn, 'w')

with open(input_fn, 'r') as f:
    line_set = set()
    for line in f:
        items = line.rstrip('\n').split('\t')
        gene = items[0].split('_')[0]
        new_line = str(gene) + '_' + str(items[1]) + '_' + str(items[2])
        line_set.add(new_line)

    for i in line_set:
        print (i, file=rmdups_outfile)
        score = float(i.split('_')[2])
        if score <= 500:
            print (i, file=scorefiltered_outfile)
