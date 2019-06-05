from functions import *
import os
import argparse
from mhc_i.src.predict_binding import Prediction

# Python 2.7


def parse_args():
    parser = argparse.ArgumentParser(description='Run mhc')
    parser.add_argument('--hla', required=True,
                        help='REQUIRED. Input the path to the hla file.')
    parser.add_argument('--patientID', required=True,
                        help='REQUIRED. Input the patient ID.')
    parser.add_argument('--output_dir', required=True,
                        help='REQUIRED. Path to output directory.')

    args = parser.parse_args()
    return args

def process_hla_mer(hla, patientID, mer_len, peptide_len, dir):
    '''
    This function processes for 1 hla and 1 mer
    hla: hla_a_03_01_01_01
    patientID: A7-A26G
    mer_len: '8'
    peptide_len: '15
    :return:
    '''


    # Make set
    IEDBList = makeSet('IEDB.txt')
    item = hla.upper().split("_")
    hla_IEDBlabel = item[0] + "-" + item[1] + "*" + item[2] + ":" + item[3]

    IEDBStr = "False"

    # new_value.append(hla_IEDBlabel)

    if hla_IEDBlabel in IEDBList:
        IEDBStr = "True"

    if IEDBStr == "False":
        hla_IEDBlabel = getClosestHLA(hla_IEDBlabel, IEDBList)

    prediction = Prediction()

    prediction.commandline_input_w_file('IEDB_recommended', hla_IEDBlabel, str(mer_len),
                                        'peptides/' + patientID + '.' + str(peptide_len) + '.txt',
                                        dir + '/IEDB_out/' + hla + '/output_IEDB.' + str(peptide_len) + '.txt')

def main():
    args = parse_args()
    hlas = []
    with open(args.hla, 'r') as f:
        for line in f:
            hlas.append(line.rstrip('\n').split('\t')[1:])
    hlas = [item for sublist in hlas for item in sublist] #make a flat list
    print hlas

    # Make output folder:

    for hla in hlas:
        os.makedirs(args.output_dir + '/IEDB_out/' + hla)

        process_hla_mer(hla, args.patientID, 8, 15, args.output_dir)
        process_hla_mer(hla, args.patientID, 9, 17, args.output_dir)
        process_hla_mer(hla, args.patientID, 10, 19, args.output_dir)
        process_hla_mer(hla, args.patientID, 11, 21, args.output_dir)

if __name__ == '__main__':
    main()

