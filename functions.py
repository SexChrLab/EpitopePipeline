import subprocess
import os
import sys
import operator
import ast
from subprocess import call
# from bs4 import BeautifulSoup
from StringIO import StringIO
import pycurl
import re
from collections import defaultdict
from collections import namedtuple
from operator import attrgetter, itemgetter
import smtplib


###################################################

def getTop200(list_):
    counter = 0
    result = []
    for val in list_:
        if counter == 2000000:
            break
        result.append(val[0])
        counter = counter + 1
    return result


###################################################

def getAnnProb(outputFilePath, num):
    key_ = []
    value_ = []
    counter = 0
    with open(outputFilePath + "Ann_output." + num + ".txt") as f:
        for line in f:
            if counter == 0:
                counter = counter + 1
                continue
            else:
                counter = counter + 1
                templine = line.split("\t")
                key_.append(templine[0].strip().replace("\"", ""))
                value_.append(float(templine[1].strip().replace("\"", "")))
    map_values = dict(zip(key_, value_))
    return map_values


###################################################

def getMapwithValues(fileName):
    key_ = []
    value_ = []
    counter = 0
    ALL_LISTS = []
    if ("output_netmhcpan.txt" in fileName):
        NetMHC = namedtuple("NetMHC", ("peptide", "start", "end", "score"))
        NetMHC_TupleList = []
        NetMHC_Set = set()
    if ("output_IEDB.txt" in fileName):
        IEDB = namedtuple("IEDB", ("peptide", "start", "end", "score"))
        IEDB_TupleList = []
        IEDB_Set = set()
    with open(fileName) as f:
        for line in f:
            if counter == 0:
                counter = counter + 1
                continue
            else:
                counter = counter + 1
                tempLine = line.split('\t')
            if ("syfpeithi.txt" in fileName):
                key_.append(tempLine[0].strip())
                value_.append(float(tempLine[1].strip()))
            if ("output_netmhcpan.txt" in fileName):
                peptideStr = tempLine[5].strip()
                if peptideStr not in NetMHC_Set:
                    NetMHC_Set.add(peptideStr)
                    key_.append(peptideStr)
                    value_.append(float(tempLine[6].strip()))
                    netmhc = NetMHC(peptide=peptideStr, start=tempLine[2].strip(), end=tempLine[3].strip(),
                                    score=tempLine[6].strip())
                    NetMHC_TupleList.append(netmhc)
            if ("output_IEDB.txt" in fileName):
                peptideStr = tempLine[5].strip()
                if peptideStr not in IEDB_Set:
                    IEDB_Set.add(peptideStr)
                    key_.append(peptideStr)
                    value_.append(float(tempLine[6].strip()))
                    iedb = IEDB(peptide=peptideStr, start=tempLine[2].strip(), end=tempLine[3].strip(),
                                score=tempLine[6].strip())
                    IEDB_TupleList.append(iedb)
    map_values = dict(zip(key_, value_))
    ALL_LISTS.append(map_values)
    if ("output_IEDB.txt" in fileName):
        ALL_LISTS.append(IEDB_TupleList)
    if ("output_netmhcpan.txt" in fileName):
        ALL_LISTS.append(NetMHC_TupleList)
    return ALL_LISTS


###################################################

def getMapwithValuesIEDB(fileName): #file name is something like output_IEDB.15.txt
    key_ = []
    value_ = []
    counter = 0
    ALL_LISTS = []

    # TODO: currently I don't see any output with netmhcpan so I'm going to comment this out to save 2 if statemetns
    # if ("output_netmhcpan" in fileName):
    #     NetMHC = namedtuple("NetMHC", ("peptide", "start", "end", "score"))
    #     NetMHC_TupleList = []
    #     NetMHC_Set = set()
    # if ("output_IEDB" in fileName):
    #     IEDB = namedtuple("IEDB", ("peptide", "start", "end", "score"))
    #     IEDB_TupleList = []
    #     IEDB_Set = set()
    IEDB = namedtuple("IEDB", ("peptide", "start", "end", "score"))
    IEDB_TupleList = []
    IEDB_Set = set()
    with open(fileName) as f:
        for line in f:
            if counter == 0:
                counter = counter + 1
                continue
            else:
                counter = counter + 1
                tempLine = line.split('\t')
            # TODO: I'm commenting out these if statemtns because right now only deal with output files with IEDB
            # if ("syfpeithi.txt" in fileName):
            #     key_.append(tempLine[0].strip())
            #     value_.append(float(tempLine[1].strip()))
            # if ("output_netmhcpan" in fileName):
            #     peptideStr = tempLine[5].strip()
            #     #  if  peptideStr not in NetMHC_Set:
            #     NetMHC_Set.add(peptideStr)
            #     key_.append(peptideStr)
            #     value_.append(float(tempLine[6].strip()))
            #     netmhc = NetMHC(peptide=peptideStr, start=tempLine[2].strip(), end=tempLine[3].strip(),
            #                     score=tempLine[6].strip())
            #     NetMHC_TupleList.append(netmhc)
            # if ("output_IEDB" in fileName):
                # peptideStr = tempLine[5].strip()
            #     # if  peptideStr not in IEDB_Set :
            #     IEDB_Set.add(peptideStr)
            #     key_.append(peptideStr)
            #     val = 0
            #     if "netmhcpan" in tempLine[6]:
            #         val = float(tempLine[14].strip())
            #         value_.append(val)
            #     elif len(tempLine) < 9:
            #         val = float(tempLine[6].strip())
            #         value_.append(val)
            #     else:
            #         val = float(tempLine[8].strip())
            #         value_.append(val)
            #     iedb = IEDB(peptide=peptideStr, start=tempLine[2].strip(), end=tempLine[3].strip(), score=val)
            #     IEDB_TupleList.append(iedb)
            peptideStr = tempLine[5].strip()
            # if  peptideStr not in IEDB_Set :
            IEDB_Set.add(peptideStr)
            key_.append(peptideStr)
            val = 0
            if "netmhcpan" in tempLine[6]:
                val = float(tempLine[14].strip())
                value_.append(val)
            elif len(tempLine) < 9:
                val = float(tempLine[6].strip())
                value_.append(val)
            else:
                val = float(tempLine[8].strip())
                value_.append(val)


            iedb = IEDB(peptide=peptideStr, start=tempLine[2].strip(), end=tempLine[3].strip(), score=val)
            IEDB_TupleList.append(iedb)
    map_values = dict(zip(key_, value_)) #key is the peptide string and valule is the score
    ALL_LISTS.append(map_values)
    # TODO: comment out these if statements
    # if ("output_IEDB" in fileName):
    #     ALL_LISTS.append(IEDB_TupleList)
    # if ("output_netmhcpan" in fileName):
    #     ALL_LISTS.append(NetMHC_TupleList)
    ALL_LISTS.append(IEDB_TupleList)
    return ALL_LISTS


###################################################

def writeNormFile(rFilePath, filepath, fileName, final_set, final_map, num):
    result = []
    key_ = []
    value_ = []
    fwrite = open(filepath + fileName + '_norm.' + num + '.csv', 'w')
    fwrite.write("Sequence," + fileName + ".bind\n")
    for val in final_set:
        fwrite.write(str(val) + "," + str(final_map[val]) + "\n")
    fwrite.close()
    rFilePath = rFilePath + "/R/"
    commandR = "Rscript " + rFilePath + fileName + "_norm.R " + filepath + " " + num
    returncode = subprocess.call(commandR, shell=True)
    counter = 0

    with open(filepath + fileName + "_normalized." + num + ".txt") as f:
        for line in f:
            if counter == 0:
                counter = counter + 1
                continue
            else:
                counter = counter + 1
                tempLine = line.split('\t')
                key_.append(tempLine[1].strip().replace("'", "").replace("\"", ""))
                value_.append(float(tempLine[3].strip()))

    map_values = dict(zip(key_, value_))

    return map_values


###################################################

def writeCombinedProb(final_output, bp_map, ann_map, ann_prob, outputFilePath, tuples, teanscriptMap_MT,
                      teanscriptMap_WT, transcript_SET):
    fwrite = open(outputFilePath + 'Epitope_prob.txt', 'w')

    for t in transcript_SET:
        value = teanscriptMap_MT[t]
        prob = bp_map[value] * (1 - ann_prob[ann_map[value]])
        peptide_tuple = [item for item in tuples if item.peptide == value]
        fwrite.write(
            "MT\t" + t + "\t" + value + "\t" + str(prob) + "\t" + peptide_tuple[0].start + "\t" + peptide_tuple[
                0].end + "\t" + peptide_tuple[0].score + "\n")
        value_WT = teanscriptMap_WT[t]
        prob = bp_map[value_WT] * (1 - ann_prob[ann_map[value_WT]])
        peptide_tuple = [item for item in tuples if item.peptide == value_WT]
        fwrite.write(
            "WT\t" + t + "\t" + value_WT + "\t" + str(prob) + "\t" + peptide_tuple[0].start + "\t" + peptide_tuple[
                0].end + "\t" + peptide_tuple[0].score + "\t" + (bp_map[value] - bp_map[value_WT]) + "\n")
    fwrite.close()


###################################################

def getANNMap(outputFilePath, num):
    result = []
    key_ = []
    value_ = []
    templine = []
    counter = 0
    fwrite = open(outputFilePath + "ANN_input_final." + num + ".txt", 'w')
    with open(outputFilePath + "Ann_input." + num + ".txt") as f:
        for line in f:
            templine = line.strip().split('\t')
            if counter == 0:
                counter = counter + 1
                if (int(num) > 14):
                    fwrite.write(
                        templine[2].strip() + "\t" + templine[3].strip() + "\t" + templine[4].strip() + "\t" + templine[
                            5].strip() + "\t" + templine[6].strip() + "\t" + templine[7].strip() + "\t" + templine[
                            8].strip() + "\t" + templine[9].strip())
                if (int(num) > 16):
                    fwrite.write("\t" + templine[10].strip())
                if (int(num) > 18):
                    fwrite.write("\t" + templine[11].strip())
                if (int(num) > 20):
                    fwrite.write("\t" + templine[12].strip())
                fwrite.write("\n")
            else:
                counter = counter + 1
                key_.append(templine[1].strip().replace("\"", ""))
                value_.append(templine[0].strip().replace("\"", ""))

                if (int(num) > 14):
                    fwrite.write(
                        templine[2].strip() + "\t" + templine[3].strip() + "\t" + templine[4].strip() + "\t" + templine[
                            5].strip() + "\t" + templine[6].strip() + "\t" + templine[7].strip() + "\t" + templine[
                            8].strip() + "\t" + templine[9].strip())
                if (int(num) > 16):
                    fwrite.write("\t" + templine[10].strip())
                if (int(num) > 18):
                    fwrite.write("\t" + templine[11].strip())
                if (int(num) > 20):
                    fwrite.write("\t" + templine[12].strip())
                fwrite.write("\n")
    fwrite.close()
    map_values = dict(zip(key_, value_))
    return map_values


###################################################

def parseSyfpeithi(motif, amers, seq):
    motif = motif.replace("H-2", "H2")
    url = "http://www.syfpeithi.de/bin/MHCServer.dll/EpitopePrediction?Motif=" + str(motif) + "&amers=" + str(
        amers) + "&SEQU=" + str(seq) + "&DoIT=++Run++"
    url = str(url)
    storage = StringIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEFUNCTION, storage.write)
    c.perform()
    c.close()
    content = storage.getvalue()
    soup = BeautifulSoup(content, 'html.parser')
    rows = soup.find_all("tr")
    # temprows = unicode(rows, "utf-8",errors="ignore")
    # rows = tables[1].find_all("tr")
    key_ = []
    val_ = []
    for row in rows[2:]:
        tds = row.find_all("td")
        counter = 0
        for td in tds:
            if counter == 0:
                counter = counter + 1
                continue
            if counter == 1:
                counter = counter + 1
                tempvar = re.sub(r'[^a-zA-Z]', "", td.getText())
                if tempvar == "gototop":
                    continue
                key_.append(tempvar.strip())
                continue
            if counter == 2:
                counter = counter + 1
                tempnum = td.getText().strip()
                val_.append(float(tempnum))
                continue
    map_values = dict(zip(key_, val_))
    return map_values


###################################################

def getSeq(file):
    seqList = []
    seq = ""
    with open(file) as f:
        for line in f:
            if ">" in line:
                continue
            else:
                seqList.append(str(line.strip()))
    seqSet = set(seqList)
    for s in seqSet:
        seq += str(s.strip())
    return seq


#############################################
# Create Sequence transcript map
#############################################
def getTranscriptSequenceMap(filepath):
    ALL_LISTS = []
    mutant_Map = defaultdict(list)
    wildType_Map = defaultdict(list)
    transcript = ""
    counter = 0
    with open(filepath) as f:
        for line in f:

            if ">" in line:
                transcript = line.strip().replace(">", "")
                trans_ = transcript.replace("MT.", "").replace("WT.", "")
                counter += 1
            else:
                seq = line.strip()
                counter += 1
            if counter % 2 == 0:

                if "MT." in transcript:
                    mutant_Map[trans_].append(seq)
                if "WT." in transcript:
                    wildType_Map[trans_].append(seq)
    ALL_LISTS.append(mutant_Map)
    ALL_LISTS.append(wildType_Map)
    return ALL_LISTS


#############################################
# Create Sequence transcript map
#############################################

def initializeDataSets(filepath):
    ALL_LISTS = []
    transcript_SET = set()
    mutant_Map = defaultdict(list)
    wildType_Map = defaultdict(list)
    transcript = ""
    counter = 0
    with open(filepath) as f:
        for line in f:
            data = line.split("\t")
            if data[0] == "MT":
                transcript_SET.add(data[1].strip())
                mutant_Map[data[1]] = data[2].strip()
            else:
                transcript_SET.add(data[1].strip())
                wildType_Map[data[1]] = data[2].strip()
    ALL_LISTS.append(transcript_SET)
    ALL_LISTS.append(mutant_Map)
    ALL_LISTS.append(wildType_Map)
    return ALL_LISTS


#############################################
# Create Sequence transcript map
#############################################

def getLowestScore(peptideList, peptideTuples):
    # print peptideList
    # print peptideTuples
    # for i in peptideTuples:
    #   #   print i
    # print len(peptideTuples)
    peptideTuples = sorted(peptideTuples, key=attrgetter('score'))
    # print peptideTuples

    for peptideTuple in peptideTuples:
        if peptideTuple.peptide in peptideList:
            return peptideTuple
            break


#############################################
# Create Sequence transcript map
#############################################  

def getFilePath(outputFilePath, hla, patientID):
    ## Some Constants
    hla_allele = hla.replace(":", "-")
    # return outputFilePath +patientID + "/" + hla_allele +  "/"
    return outputFilePath + "/" + hla_allele + "/"


###################################################

def makeSet(filename):
    with open(filename, 'r') as f:
        return set([line.strip() for line in f])


# def readSyfpeithi():
#     with open("syfpeithi.txt") as f: #TODO: Give the full path to this file?
#         syfpeithi = [line.strip() for line in f]
#     return syfpeithi
#
# def readIEDB():
#      IEDB = []
#      with open("IEDB.txt" ) as f: #TODO: Give the full path to this file?
#           for line in f:
#                IEDB.append(line.strip())
#      return IEDB
#
# def readnetmhcpan():
#      netmhcpan = []
#      with open("netmhcpan.txt" ) as f:
#           for line in f:
#                netmhcpan.append(line.strip())
#      return netmhcpan

def writeInputFile(filepath, file, patientID, num):
    seq = ""
    header = ""
    data = ""
    filename = filepath + patientID + "." + num + '.txt'
    # peptideFileName = filepath  + file
    peptideFileName = filepath + num + "mers/" + file
    fwrite = open(filename, 'w')
    with open(peptideFileName) as f:
        for line in f:
            if ">" in line:
                header = line.strip();
                data = ""
            else:
                data = line.strip().replace("X", "")
            n = int(num) * 2
            if len(data) > (int(num) / 2):
                # print header + "\n" + data+"\n"
                fwrite.write(header + "\n" + data + "\n")

    fwrite.close();


def getClosestHLA(hla, strList):
    minDiff = 1000000
    minDiffHLa = ""
    hlaLists = [s for s in strList if hla[:6] in s]
    for h in hlaLists:
        diff = abs(float(hla[6:].replace(":", "")) - float(h[6:].replace(":", "")))
        if (diff <= minDiff):
            minDiff = diff
            minDiffHLa = h
    return minDiffHLa


def getAllSeq(seq, length): #This function takes in a sequence such as CWVRDRNLRPKFSQI and length of 8 #TODO: refactor this function
    l = len(seq) + 1 #for example l = 15+1 = 16
    list_ = []
    for i in range(l - length): #range(8)
        list_.append(seq[: length])
        seq = seq[1:]
    return list_


#############################################
# Create Sequence transcript map
#############################################

def getMutantWildTypeData(filepath, patientId,
                          num):  # this function takes in the filepath, patientID and either 15, 17, 19, or 21
    ALL_LISTS = []
    # seq_transcript_map = defaultdict(list) #this dict doesn't get returned anywhere
    transcript_SET = set()

    filename = patientId + "." + str(num) + ".txt"
    transcript = ""
    trans_ = ""
    counter = 0
    mutant_Map = defaultdict(list) #Key is the transcript name and values are all the possible x-mers
    wildType_Map = defaultdict(list) #Key is the transcript name and values are all the possible x-mers
    with open(filepath + "peptides/" + filename) as f:  # open A7-A26G.15.txt
        for line in f:
            seq = ""
            if ">" in line:
                transcript = line.strip().replace(">", "")
                trans_ = transcript.replace("MT.", "").replace("WT.", "")
                transcript_SET.add(trans_)
                counter += 1
            else:
                seq = line.strip()
                counter += 1
            if counter % 2 ==0:
                length = (int(num) / 2) + 1
                # seq_transcript_map[seq].append(transcript)
                if "MT." in transcript:
                    mutant_Map[trans_].append(getAllSeq(seq, length))
                    # s = mutant_Map[trans_]
                    # fwriteMT.write("MT\t" + str(trans_) + "\t" + str(s[0]) + "\n")
                if "WT." in transcript:
                    wildType_Map[trans_].append(getAllSeq(seq, length))

    ALL_LISTS.append(transcript_SET)
    ALL_LISTS.append(mutant_Map)
    ALL_LISTS.append(wildType_Map)
    return ALL_LISTS


def getSameSeqScore(mutantLowestScoreTuple, wildTypePeptideList, peptideTuples):
    peptide_tuple_sameSeq = [item for item in peptideTuples if str(item.start) == str(mutantLowestScoreTuple.start)]

    for tup in peptide_tuple_sameSeq:
        if tup.peptide in wildTypePeptideList:
            return tup


def getPeptides(IEDB_transcriptMap, num):
    peptideSet = set()
    for item in IEDB_transcriptMap:
        if item.mer == num:
            if item.data == None:
                peptideSet.add(item.data.peptide)
    return peptideSet


def writeTransformFile(filepath, outputFilePath, norm_map, peptideSet, num):
    bp_key_ = []
    bp_value_ = []
    fwrite = open(outputFilePath + 'transform_input.' + num + '.csv', 'w')
    fwrite.write("Epitope,IEDB.Norm,BP Score\n")
    final_set = set(peptideSet)
    for val in final_set:
        val = val.replace("'", "")
        bp_score = norm_map[val];
        bp_key_.append(val)
        bp_value_.append(bp_score)
        fwrite.write(val + "," + str(norm_map[val]) + "," + str(bp_score) + "\n")

    fwrite.close()
    bp_map = dict(zip(bp_key_, bp_value_))
    rFilePath = filepath + "/R/"
    commandR = "Rscript " + rFilePath + "/Peptidematrix.R " + outputFilePath + " " + rFilePath + " " + num
    returncode = subprocess.call(commandR, shell=True)
    return bp_map


def getData(dataMapArr, num):
    dataTuples = namedtuple("dataTuple", ("mer", "data"))
    dataTuples = dataMapArr
    t = list()
    for s in dataTuples:
        if s.mer == num:
            dataField = namedtuple("IEDB", ("peptide", "start", "end", "score"))
            dataField = s.data
            t.append(str(dataField.start))
            t.append(str(dataField.peptide))
            t.append(str(dataField.score))
    return tuple(t)


def getNewHLA(newOldHLAMap, hla):
    if hla == newOldHLAMap[hla]:
        return ""
    else:
        isSameHla = newOldHLAMap[hla]
        isSameHla = isSameHla.replace(":", "-")
        return isSameHla


def getFinalMer(merList):
    mer = sorted(merList, key=itemgetter(1))
    return mer[0]


def readSyfHLAfile(filepath):
    key_ = []
    value_ = []

    with open(filepath + "/hla-syf.txt") as f:
        for line in f:
            line = line.split('\t')
            key_.append(line[0])
            value_.append(line[1].strip())
    return dict(zip(key_, value_))


def getSeqScores(peptideList, peptideTuples):
    key_ = []
    value_ = []
    for peptideTuple in peptideTuples:
        if peptideTuple.peptide in peptideList:
            key_.append(peptideTuple.peptide)
            value_.append(peptideTuple.score)
    return dict(zip(key_, value_))


def getlowestScore3(peptideList, IEDB_map, NetMHC_map, Syfp_map, t):
    val_ = []
    key_ = []
    tab = "\t"
    for peptide in peptideList:
        x1 = [v for i, v in enumerate(IEDB_map) if v[0] == peptide]
        peptide, IEDB_score = x1[0]

        x2 = [v for i, v in enumerate(NetMHC_map) if v[0] == peptide]
        peptide, NetMhc_score = x2[0]
        if len(Syfp_map) > 0:
            x3 = [v for i, v in enumerate(Syfp_map) if v[0] == peptide]
            peptide, Syfp_score = x3[0]

            val_.append((IEDB_score + NetMhc_score + Syfp_score) / 3)
        else:
            val_.append((IEDB_score + NetMhc_score) / 2)
        key_.append(peptide)

    map_values = dict(zip(key_, val_))
    # print map_values
    lowest_value = [sorted(map_values.items(), key=lambda x: x[1])][0]
    lowest_value = lowest_value[0]
    # print lowest_value
    return lowest_value


def getDataIEDB(dataMapArr, num, transcript):
    # print dataMapArr
    dataTuples = namedtuple("dataTuple", ("mer", "data"))
    dataTuples = dataMapArr
    t = list()
    for s in dataTuples:
        if s.mer == num:
            # dataField = namedtuple("IEDB", ("peptide","start","end","score"))
            dataField = s.data
            t.append(str(dataField[0]))
            t.append(str(dataField[1]))
            # t.append(str(dataField.score))
    return tuple(t)


def getHeaderText():
    return "PatientId" + '\t' + "Allele" + '\t' + "Gene-transcript" + '\t' + \
           "Mutant" + '\t' + "Start_15mer" + '\t' + "Peptide_15mer" + '\t' + "IEDB_Binding_15mer" + '\t' + \
           "Start_Position_17mer" + '\t' + "Peptide_17mer" + '\t' + "IEDB_Binding_17mer" + '\t' + \
           "Start_Position_19mer" + '\t' + "Peptide_19mer" + '\t' + "IEDB_Binding_19mer" + '\t' + \
           "Start_Position_21mer" + '\t' + "Peptide_21mer" + '\t' + "IEDB_Binding_21mer" + '\t' + \
           "BestPeptide_Mutant" + '\t' + "BestScore_Mutant" + '\t' + "BestScoreMer_Mutant" + '\t' + \
           "WildType" + '\t' + "Start_15mer" + '\t' + "Peptide_15mer" + '\t' + "IEDB_Binding_15mer" + '\t' + \
           "Start_Position_17mer" + '\t' + "Peptide_17mer" + '\t' + "IEDB_Binding_17mer" + '\t' + \
           "Start_Position_19mer" + '\t' + "Peptide_19mer" + '\t' + "IEDB_Binding_19mer" + '\t' + \
           "Start_Position_21mer" + '\t' + "Peptide_21mer" + '\t' + "IEDB_Binding_21mer" + '\t' + \
           "BestPeptide_WildType" + '\t' + "BestScore_WildType" + '\t' + "BestScoreMer_WildType" + '\t' + \
           "Same_Seq" + '\t' + "Start_15mer" + '\t' + "Peptide_15mer" + '\t' + "IEDB_Binding_15mer" + '\t' + \
           "Start_Position_17mer" + '\t' + "Peptide_17mer" + '\t' + "IEDB_Binding_17mer" + '\t' + \
           "Start_Position_19mer" + '\t' + "Peptide_19mer" + '\t' + "IEDB_Binding_19mer" + '\t' + \
           "Start_Position_21mer" + '\t' + "Peptide_21mer" + '\t' + "IEDB_Binding_21mer" + '\t' + \
           "BestPeptide_SameSeq" + '\t' + "BestScore_SameSeq" + '\t' + "BestScoreMer_SameSeq" + '\t' + \
           "isSameHLA\n"


def getActualScores(peptide, IEDB, NetMhc, Syfp):
    data = []
    if len(peptide) == 8:
        if peptide in IEDB[0]:
            data.append(IEDB[0][peptide])
        if peptide in NetMhc[0]:
            data.append(NetMhc[0][peptide])
        if Syfp and peptide in Syfp[0]:
            data.append(Syfp[0][peptide])

    if len(peptide) == 9:
        if peptide in IEDB[1]:
            data.append(IEDB[1][peptide])
        if peptide in NetMhc[1]:
            data.append(NetMhc[1][peptide])
        if Syfp and peptide in Syfp[1]:
            data.append(Syfp[1][peptide])

    if len(peptide) == 10:
        if peptide in IEDB[2]:
            data.append(IEDB[2][peptide])
        if peptide in NetMhc[2]:
            data.append(NetMhc[2][peptide])
        if Syfp and peptide in Syfp[2]:
            data.append(Syfp[2][peptide])

    if len(peptide) == 11:
        if peptide in IEDB[3]:
            data.append(IEDB[3][peptide])
        if peptide in NetMhc[3]:
            data.append(NetMhc[3][peptide])
        if Syfp and peptide in Syfp[3]:
            data.append(Syfp[3][peptide])

    return data
