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
                        help='REQUIRED. Input the patient ID. PatientID is also the name of the output directory')
    parser.add_argument('--hla_type', required=True,
                        help='REQUIRED.')

    args = parser.parse_args()
    return args

def process_hla_mer(hla, patientID, mer_len, peptide_len, hla_type):
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
                                        patientID + '/peptides/' + patientID + '.' + str(peptide_len) + '.txt',
                                        patientID + '/IEDB_out/' + hla_type + '/output_IEDB.' + str(peptide_len) + '.txt')

def main():
    args = parse_args()
    # hlas = []

    # The following codes commented out are for the hla file format from polysolver
    # with open(args.hla, 'r') as f:
    #     for line in f:
    #         hlas.append(line.rstrip('\n').split('\t')[1:])
    # hlas = [item for sublist in hlas for item in sublist] #make a flat list
    # print hlas

    # The following codes are for the hla file format from HLA-LA (after running the script format_hla_output.py file)
    with open(args.hla, 'r') as f:
        for line in f:
            hla = line.rstrip('\n').split('\t')[1]
            # hlas.append(line.rstrip('\n').split('\t')[1])
            # hla = line.rstrip('\n').split('\t')[1]
            print hla

    # Make output folder:

    # for hla in hlas:
    #     if not os.path.exists(args.output_dir + '/IEDB_out/' + hla):
    #         os.makedirs(args.output_dir + '/IEDB_out/' + hla)

            if not os.path.exists(args.patientID + '/IEDB_out/' + args.hla_type):
                os.makedirs(args.patientID + '/IEDB_out/' + args.hla_type)

            process_hla_mer(hla, args.patientID, 8, 15, args.hla_type)
            process_hla_mer(hla, args.patientID, 9, 17, args.hla_type)
            process_hla_mer(hla, args.patientID, 10, 19, args.hla_type)
            process_hla_mer(hla, args.patientID, 11, 21, args.hla_type)

if __name__ == '__main__':
    main()

