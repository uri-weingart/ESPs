import sys
import json
import Bio
import datetime
import platform
import string
import argparse
import treelib
from treelib import Node, Tree
import pickle
import ahocorasick
import Bio
from Bio import SeqIO
##########################################################################
#   Decode input parms                                                                        #
##########################################################################


def read_params(x):
    CommonArea = dict()
    parser = argparse.ArgumentParser(
        description='Search Specific Peptides (SPs)  within proteins sequences')

    parser.add_argument(
        '-i',
        '--InputFasta',
        action="store",
        dest='InputFasta',
        nargs='?',
        help='Input proteins file name (Format: fasta)   Default = Input.fasta',
        required=False,
        default="data/PTest.fasta")

    parser.add_argument(
        '-dSPs',
        '--Input_dSPs_json',
        action="store",
        dest='dSPs',
        nargs='?',
        help='json file containing all the SPs (Provided by the package) - Default =  ESPs_Augmented.json',
        required=False,
        default="ESPs_Augmented.json")

    parser.add_argument(
        '-Au',
        '--Automaton_pickle',
        action="store",
        dest='Automaton_pickle',
        nargs='?',
        help='Pickle automaton required for aho_corasick search of multiple SPs in a sequence',
        required=False,
        default="ESPs_Automaton_Pickle_Augmented")


    parser.add_argument('-o', '--Output',
                        action="store",
                        dest='Output_FileName',
                        nargs='?',
                        help='Output file name',
                        required=False,
                        default="Hits_Summary.csv")




    parser.add_argument(
        '--SP_Len_Thresh',
        action="store",
        dest='SP_Len_Thresh',
        help='Threshold minumum coverage for a prediction - Default is 7 ** Very important parameter ** ',
        nargs='?',
        required=False,
        default="7")

    CommonArea['parser'] = parser
    return CommonArea


##########################################################################
#   Prepare run parameters                                                                    #
##########################################################################
def Prepare_Run_Parameters(results):
    dParms = dict()
    dParms["Input_Fasta_Name"] = results.InputFasta
    dParms["dSPs_File_Name"] = results.dSPs
    dParms["Output_FileName"] = results.Output_FileName
    dParms["Automaton_Pickle_File"] = results.Automaton_pickle

    Check_Number = results.SP_Len_Thresh
    if Check_Number.isnumeric():
        if 5 <= int(Check_Number) <= 100:
            dParms["SP_Len_Thresh"] = int(Check_Number)

    sDashes = "-" * 60
    print(sDashes)
    print("*  Parameters chosen: ")
    for sParm in sorted(dParms.keys()):
        sOut = "*  " + sParm + ": " + str(dParms[sParm]) + " "
        print(sOut)

    print(sDashes)
    return (dParms)



#******************************************************************************
#     Predict EC from Tree                                                    *
#******************************************************************************
def Predict_EC_From_Tree(tree):
 
	lNodes_Level1 = tree.children("Root")
	lEC_Predictions = list()
	for node in lNodes_Level1:
		sub_tree = tree.subtree(node.identifier)
		sub_tree_Depth = sub_tree.depth()
		lSub_Tree  =  [sub_tree[node].tag for node in  sub_tree.expand_tree(mode=Tree.WIDTH)]
		lECs = [sub_tree.level(node) for node in lSub_Tree]
		for Level in range(sub_tree.depth() + 1):
			Entries = sub_tree.size(level = Level)
			if Entries > 1:
				iPrediction_Level = Level - 1
				iEC_Prediction = lECs.index(iPrediction_Level)
				EC_Prediction = lSub_Tree[iEC_Prediction]
				lEC_Predictions.append(EC_Prediction)
				break
			if Entries  == 1 and Level == sub_tree_Depth:
				iPrediction_Level = Level 
				iEC_Prediction = lECs.index(iPrediction_Level)
				EC_Prediction = lSub_Tree[iEC_Prediction]
				lEC_Predictions.append(EC_Prediction)
				break
	return(lEC_Predictions)

# ******************************************************************************
#     Add  EC to the tree                                                     *
# ******************************************************************************


def Add_EC_To_Tree(tree, EC,  Coverage):


    EC = EC.replace(".-", "")
    EC_Level = EC.count(".") + 1

    if EC_Level == 1:
        ECLevel1 = EC.split(".")[0]
        if not tree.contains(ECLevel1):
            tree.create_node(ECLevel1, ECLevel1, parent="Root", data=Coverage)
    if EC_Level == 2:
        ECLevel1 = EC.split(".")[0]
        ECLevel2 = EC.split(".")[0] + "." + EC.split(".")[1]
        if not tree.contains(ECLevel1):
            tree.create_node(ECLevel1, ECLevel1, parent="Root", data=0)
        tree.create_node(ECLevel2, ECLevel2, parent=ECLevel1, data=Coverage)

    if EC_Level == 3:
        ECLevel1 = EC.split(".")[0]
        ECLevel2 = EC.split(".")[0] + "." + EC.split(".")[1]
        ECLevel3 = EC.split(".")[0] + "." + \
            EC.split(".")[1] + "." + EC.split(".")[2]
        if not tree.contains(ECLevel1):
            tree.create_node(ECLevel1, ECLevel1, parent="Root", data=0)
        if not tree.contains(ECLevel2):
            tree.create_node(ECLevel2, ECLevel2, parent=ECLevel1, data=0)
        tree.create_node(ECLevel3, ECLevel3, parent=ECLevel2, data=Coverage)

    if EC_Level == 4:
        ECLevel1 = EC.split(".")[0]
        ECLevel2 = EC.split(".")[0] + "." + EC.split(".")[1]
        ECLevel3 = EC.split(".")[0] + "." + \
            EC.split(".")[1] + "." + EC.split(".")[2]
        ECLevel4 = EC.split(".")[0] + "." + EC.split(".")[1] + \
            "." + EC.split(".")[2] + "." + EC.split(".")[3]
        if not tree.contains(ECLevel1):
            tree.create_node(ECLevel1, ECLevel1, parent="Root", data=0)
        if not tree.contains(ECLevel2):
            tree.create_node(ECLevel2, ECLevel2, parent=ECLevel1, data=0)
        if not tree.contains(ECLevel3):
            tree.create_node(ECLevel3, ECLevel3, parent=ECLevel2, data=0)
        if not tree.contains(ECLevel4):
            tree.create_node(
                ECLevel4,
                ECLevel4,
                parent=ECLevel3,
                data=Coverage)


    return (tree)


# ******************************************************************************
#     Calculate Hits coverage                                                 *
# ******************************************************************************
def Map_SP_Hits_by_EC(lSP_Hits, dSPs):
    dCoverage = dict()
    dHits_by_EC = dict()
    dMap = dict()
    for ent in lSP_Hits:
        SP = ent[0]
        EC = dSPs[SP]['EC']
        lHits = ent[1]
        if EC not in dHits_by_EC:
            dHits_by_EC[EC] = list()
        for i_hit_location in lHits:
            lHit_Entry = [i_hit_location, i_hit_location + len(SP) - 1]
            dHits_by_EC[EC].append(lHit_Entry)
    for EC in dHits_by_EC.keys():
        dMap[EC] = dict()
        for ent in dHits_by_EC[EC]:
            iFrom = ent[0]
            i_To = ent[1]
            for indx in range(iFrom, i_To + 1):
                if indx not in dMap[EC]:
                    dMap[EC][indx] = 0
                dMap[EC][indx] += 1
    for EC in sorted(dMap.keys()):
        Coverage = 0
        for AA_loc in dMap[EC].keys():
            Coverage+=dMap[EC][AA_loc]
        dCoverage[EC] = Coverage
    return dCoverage


#******************************************************************************
#    Accumulate coverage up in the tree (From leaves up)                      *
#******************************************************************************
def    Accumulate_up(tree):


    lTree = tree.expand_tree(mode=1)
 

    lECs = [[],[],[],[],[]]

    for identifier in lTree:
        depth = tree.depth(identifier) 
        tree_entry = tree.get_node(identifier)
        lECs[depth].append(tree_entry)


    for EC_Lvl in (3,2,1,0):
        for entry in lECs[EC_Lvl]:
            Accum_Coverage = 0

            lChildren = tree.children(entry.identifier)
            for child in lChildren:
                Accum_Coverage+=child.data

            tree.update_node(entry.identifier, data = Accum_Coverage + entry.data)
                
   

    return tree 

# ******************************************************************************
#     Program Starts here                                                      *
# ******************************************************************************
StartTime = datetime.datetime.now()
CommonArea = read_params(sys.argv)
parser = CommonArea['parser']
results = parser.parse_args()
dParms = Prepare_Run_Parameters(results)
cT ="\t"
cN = "\n"




# ******************************************************************************
#     Load the SP trie and build automaton                                    *
# ******************************************************************************
pickle_filename = dParms["Automaton_Pickle_File"]
infile = open(pickle_filename, 'rb')
A_SPs = pickle.load(infile)
infile.close()
A_SPs.make_automaton()

End_Load_Automaton_and_SPs_Time = datetime.datetime.now()



with open(dParms["dSPs_File_Name"]) as json_file:
    dSPs = json.load(json_file)
Input_Fasta_File = dParms["Input_Fasta_Name"]

outfile = open(dParms["Output_FileName"], "w")

myfasta_sequences = SeqIO.parse(open(Input_Fasta_File), 'fasta')
for seq in myfasta_sequences:
    dNew_Hits = dict()
    lNew_Hits = list()
    Seq = str(seq.seq)
    iLenSeq = len(Seq)
    AC = seq.id

    tree = None
    tree = Tree()
    tree.create_node("Root", "Root", data=0)  # root node



    for end_index, (insert_order, original_value) in A_SPs.iter(Seq):
        start_index = end_index - len(original_value) + 1
        if original_value not in dNew_Hits:
            dNew_Hits[original_value] = list()
        dNew_Hits[original_value].append(start_index + 1)

    lSP_Hits = list()
    for original_value in sorted(dNew_Hits):
        lNew_Entry = [original_value]
        lNew_Entry.append(dNew_Hits[original_value])
        lSP_Hits.append(lNew_Entry)

    dCoverage = Map_SP_Hits_by_EC(lSP_Hits, dSPs)



    for EC in sorted(dCoverage.keys()):
        Coverage  =  dCoverage[EC]
        tree = Add_EC_To_Tree(tree, EC,   Coverage) 

 


    tree  =  Accumulate_up(tree)  # Accumulate the coverage up in the tree, from leaves up the branches 



    lEC_Predictions  = Predict_EC_From_Tree(tree) 
    Outrec = AC + cT 
    if len(lEC_Predictions) == 0:
        strMgs = " No Prediction  ==> "  
    else:
        strMgs = " Predictions: ==>  " +      '~'.join(EC for EC in lEC_Predictions)  
    Outrec = Outrec + strMgs + cN
    outfile.write(Outrec)


outfile.close()

EndTime = datetime.datetime.now()
print("Start Time                               : ", StartTime.strftime("%Y-%m-%d %H:%M:%S"))
print("Ended loading tables and Start Searching : ", End_Load_Automaton_and_SPs_Time.strftime("%Y-%m-%d %H:%M:%S"))
print("End Time                                 : ", EndTime.strftime("%Y-%m-%d %H:%M:%S"))
print("Program Completed Successfully")
sys.exit(0)
