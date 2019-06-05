# Run script `` in order to run netMHC to produce file such as output_IEDB.15.txt.
# This script is intended to process after IEDB step.
# It will:
# 1. Get min score
# 2. Combine 15mer, 17mer, 19mer, and 21mer
# 3. Find the mer with the lowest IEBD score
# 4. Find potential neoepitope (IEDB <= 500)

from collections import defaultdict
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Find potential neoepitope')
    parser.add_argument('--dir', required=True,
                        help='REQUIRED. Input the path to the directory.')
    parser.add_argument('--hla', required=True,
                        help='REQUIRED. Input the path to the hla file.')
    parser.add_argument('--patientID', required=True,
                        help='REQUIRED. Input the patient ID.')
    args = parser.parse_args()
    return args

def get_min_score(mer_len, output_IEDB_fn, peptide_fn, min_score_out_fn):
    score_dict = defaultdict(list)  # key is the peptide and values are scores
    with open(output_IEDB_fn, 'r') as f:
        for line in f:
            if not line.startswith('allele'):
                tempLine = line.rstrip('\n').split('\t')
                if "netmhcpan" in tempLine[6]:
                    score_dict[tempLine[5]].append(float(tempLine[14]))
                elif len(tempLine) < 9:
                    score_dict[tempLine[5]].append(float(tempLine[6]))
                else:
                    score_dict[tempLine[5]].append(float(tempLine[8]))

    transcript = ""
    trans_ = ""
    counter = 0
    outfile = open(min_score_out_fn, 'w')
    with open(peptide_fn, 'r') as f:  # open A7-A26G.15.txt
        for line in f:
            seq = ""
            if ">" in line:
                transcript = line.strip().replace(">", "")
                trans_ = transcript.replace("MT.", "").replace("WT.", "")
                counter += 1
            else:
                seq = line.strip()
                counter += 1
            if counter % 2 == 0:
                if "MT." in transcript:
                    peptideScore = defaultdict()
                    all_mer_len_peptide = [seq[i:i + int(mer_len)] for i in range(int(mer_len))]
                    for mer_len_peptide in all_mer_len_peptide:
                        if len(score_dict[mer_len_peptide]) != 0:
                            peptideScore[mer_len_peptide] = min(score_dict[mer_len_peptide])
                    # print peptideScore
                    out = min(peptideScore.items(), key=lambda l: l[1])
                    out_to_print = [trans_, out[0], str(out[1])]
                    print ('\t'.join(x for x in out_to_print), file=outfile)


def select_peptide(merged_peptides_fn, out_fn):
    outfile = open(out_fn, 'w')
    with open(merged_peptides_fn, 'r') as f:
        for line in f:
            if not line.startswith('transcript'):
                peptides_scores = {}
                items = line.rstrip('\n').split('\t')
                for i in (1, 3, 5, 7):
                    peptides_scores[items[i]] = float(items[i+1])
                lowest_peptide_score = min(peptides_scores.items(), key=lambda l:l[1])
                out = [items[0], lowest_peptide_score[0], str(lowest_peptide_score[1])]
                print ('\t'.join(x for x in out), file=outfile)

def main():
    args = parse_args()

    hlas = []
    with open(args.hla, 'r') as f:
        for line in f:
            hlas.append(line.rstrip('\n').split('\t')[1:])
    hlas = [item for sublist in hlas for item in sublist] #make a flat list

    for hla in hlas:

        get_min_score(8, args.dir + '/IEDB_out/' + hla + '/output_IEDB.15.txt', args.dir + '/peptides/'+ args.patientID + '.15.txt', args.dir + '/IEDB_out/' + hla + '/output_IEDB.15_minscore.txt')
        get_min_score(9, args.dir + '/IEDB_out/' + hla + '/output_IEDB.17.txt',
                      args.dir + '/peptides/' + args.patientID + '.17.txt',
                      args.dir + '/IEDB_out/' + hla + '/output_IEDB.17_minscore.txt')
        get_min_score(10, args.dir + '/IEDB_out/' + hla + '/output_IEDB.19.txt',
                      args.dir + '/peptides/' + args.patientID + '.19.txt',
                      args.dir + '/IEDB_out/' + hla + '/output_IEDB.19_minscore.txt')
        get_min_score(11, args.dir + '/IEDB_out/' + hla + '/output_IEDB.21.txt',
                      args.dir + '/peptides/' + args.patientID + '.21.txt',
                      args.dir + '/IEDB_out/' + hla + '/output_IEDB.21_minscore.txt')

        peptide_15_df = pd.read_csv(args.dir + '/IEDB_out/' + hla + '/output_IEDB.15_minscore.txt', sep='\t', names = ["transcript", "peptide_15", "score_15"])
        peptide_17_df = pd.read_csv(args.dir + '/IEDB_out/' + hla + '/output_IEDB.17_minscore.txt', sep='\t',
                                    names=["transcript", "peptide_17", "score_17"])
        peptide_19_df = pd.read_csv(args.dir + '/IEDB_out/' + hla + '/output_IEDB.19_minscore.txt', sep='\t',
                                    names=["transcript", "peptide_19", "score_19"])
        peptide_21_df = pd.read_csv(args.dir + '/IEDB_out/' + hla + '/output_IEDB.21_minscore.txt', sep='\t',
                                    names=["transcript", "peptide_21", "score_21"])

        peptide_15_17_merge = pd.merge(peptide_15_df, peptide_17_df, on="transcript")
        peptide_15_17_19_merge = pd.merge(peptide_15_17_merge, peptide_19_df, on="transcript")
        peptide_15_17_19_21_merge = pd.merge(peptide_15_17_19_merge, peptide_21_df, on="transcript")

        peptide_15_17_19_21_merge.to_csv(args.dir + '/IEDB_out/' + hla + '/output_IEDB.15.17.19.21_minscore.txt', sep='\t', index=False)

        select_peptide(args.dir + '/IEDB_out/' + hla + '/output_IEDB.15.17.19.21_minscore.txt', args.dir + '/IEDB_out/' + hla + '/output_IEDB.15.17.19.21_minscore_minpeptide.txt')

main()


# # arguments should be the length of the mer and hla types
# mer_len = int(sys.argv[1]) #for example 8
# peptide_len = int(sys.argv[2]) #for example 15

# Obtain a dictionary of score





