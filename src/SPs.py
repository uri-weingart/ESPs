import sys
import pdb
import json
import Bio
from Bio import SeqIO
import datetime
import platform
import re
import pandas as pd
import numpy as np
import string
import collections
from collections import Counter
import argparse
from toolz import unique
import operator
from functools import reduce
import treelib
from treelib import Node, Tree
 
###############################################################################################
#   Decode input parms                                                                        #
###############################################################################################
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Search Specific Peptides (SPs)  within proteins sequences')
	
	parser.add_argument('-i',   '--InputFasta' ,  
	   action = "store",
	   dest = 'InputFasta',
	   nargs = '?',
	   help='Input proteins file name (Format: fasta)   Default = Input.fasta',
	   required = False,
	   default="data/test.fasta")
	   
	parser.add_argument('-dSPs',   '--Input_dSPs_json' ,  
	   action = "store",
	   dest = 'dSPs',
	   nargs = '?',
	   help='json file containing all the SPs (Provided by the package) - Default =  dSPs.json',
	   required = False,
	   default = "ESPs.json")

	parser.add_argument('-o',   '--Output' ,  
	   action = "store",
	   dest = 'Output_Html_FileName',
	   nargs = '?',
	   help='Output html file: Default = Hits_Summary.html',
	   required = False,
	   default = "Hits_Summary.html")

	parser.add_argument('--SP_Len_Thresh',
	   action="store",
	   dest='SP_Len_Thresh',
	   help='Threshold minumum coverage for a prediction - Default is 7 ** Very important parameter ** ', 
	   nargs='?',
	   required = False,
	   default = "7")	  

	parser.add_argument('--Print_Html_Report' ,  
	   action = "store",
	   dest = 'Print_Html_Report',
	   nargs = '?',
	   help='Flag: Y means print the main html report',
	   required = False,
	   default = "Y")	

	parser.add_argument('--w_Seq_Display',
	   action="store",
	   dest='w_Seq_Display',
	   help='Width of line containing sequence display- Default = 80  Valid: between 30 and 100', 
	   nargs='?',
	   required = False,
	   default = "80")

	parser.add_argument('--Output_xlsx' ,  
	   action = "store",
	   dest = 'Output_xlsx_FileName',
	   nargs = '?',
	   help='Output Excel file (Optional): Default = Hits_Summary.xlsx',
	   required = False,
	   default = "Hits_Summary.xlsx")		   
	
	parser.add_argument('--Output_xlsx_Required' ,  
	   action = "store",
	   dest = 'Output_xlsx_Required',
	   nargs = '?',
	   help='Flag: Y means Create excel output file -  works in conjunction with --Output_xlsx,  Default = N ',
	   required = False,
	   default = "N")	
	   
	parser.add_argument('--Output_Simple_Html' ,  
	   action = "store",
	   dest = 'Output_Simple_Html_FileName',
	   nargs = '?',
	   help='Output in simple html (just the hits analysis and no colors) - Default:  Hits_Simple.html',
	   required = False,
	   default = "Hits_Simple.html")		   
	
	parser.add_argument('--Output_Simple_Html_Required' ,  
	   action = "store",
	   dest = 'Output_Simple_Html_Required',
	   nargs = '?',
	   help='Flag: Y means Create simple html  output file -  woorks in conjunction with ----Output_Simple_Html,  Default = N ',
	   required = False,
	   default = "N")	
	   
	parser.add_argument('--Output_Predictions_xlsx_Required' ,  
	   action = "store",
	   dest = 'Output_Predictions_xlsx_Required',
	   nargs = '?',
	   help='Flag: Y means Create excel output file contaning only the predictions -  woorks in conjunction with --Output_xlsx,  Default = N ',
	   required = False,
	   default = "N")	   

	parser.add_argument('--Output_Predictions_xlsx' ,  
	   action = "store",
	   dest = 'Output_Predictions_xlsx_FileName',
	   nargs = '?',
	   help='Output Excel file containing the predictions (Optional): Default = Summary_Predictions.xlsx',
	   required = False,
	   default = "Summary_Predictions.xlsx")	   
	
	CommonArea['parser'] = parser
	return  CommonArea


###############################################################################################
#   Prepare run parameters                                                                    #
###############################################################################################
def Prepare_Run_Parameters(results):
		dParms = dict()
		dParms["Input_Fastq_Name"]  = results.InputFasta
		dParms["dSPs_File_Name"]  = results.dSPs	
		dParms["Output_Html_FileName"]  = results.Output_Html_FileName

		Check_Number = results.SP_Len_Thresh
		if Check_Number.isnumeric():
			if  5 <=  int(Check_Number)  <= 100:
					dParms["SP_Len_Thresh"] =  int(Check_Number )
		
		dParms["w_Seq_Display"] = 80
		Check_Number = results.w_Seq_Display
		if Check_Number.isnumeric():
			if  30 <=  int(Check_Number)  <= 100:
					dParms["w_Seq_Display"] =  int(Check_Number )
				
		dParms["Output_xlsx_Required"] = results.Output_xlsx_Required
		dParms["Output_xlsx_FileName"] = results.Output_xlsx_FileName
		
		dParms["Output_Pred_xlsx_Required"] = results.Output_Predictions_xlsx_Required
		dParms["Output_Predictions_xlsx_FileName"] = results.Output_Predictions_xlsx_FileName
		
		dParms["Output_Simple_Html_Required"] = results.Output_Simple_Html_Required
		dParms["Output_Simple_Html_FileName"] = results.Output_Simple_Html_FileName
		dParms["Print_Html_Report"] = results.Print_Html_Report

		sDashes = "-"*60
		print(sDashes)
		print("*  Parameters chosen: ")
		for sParm in  sorted(dParms.keys()):
			sOut = "*  " + sParm + ": " + str(dParms[sParm]) + " "
			print(sOut)
		
		print(sDashes)
		return( dParms )

#******************************************************************************
#     Display coverage                                                        *
#******************************************************************************
def  Display_Coverage(AC, dACs_Hits, FileOut_Html):
	lDisplay_Output = list()
	rc = Write_Rec_to_Html(4, '<table> ' , FileOut_Html)
	for ind  in range( len( dACs_Hits[AC]["lProperties"])):
		dProperty = dACs_Hits[AC]["lProperties"][ind]
		for Property in dProperty.keys():
			Property_Value = dProperty[Property]
			Property_Coverage = dACs_Hits[AC]["Coverage"][ind]
			lDisplay_Entry = [Property, Property_Value, Property_Coverage]
			lDisplay_Output.append(lDisplay_Entry)
		lDisplay_Output1 = sorted(lDisplay_Output)
	
	sHeader_Text = "Coverage analysis of " + dACs_Hits[AC] ["Description"]
	sHeader = '<h3 class="WhiteText">' + 	sHeader_Text + "</h3>"
	Write_Rec_to_Html(0,sHeader , FileOut_Html)
	rc = Write_Rec_to_Html(0, '<table border="1">' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '<tr style="text-align: center;">' , FileOut_Html)
	lHeader = ["Function", "Function Value", "Coverage"]
	for ent in lHeader:
		sOut = "<th>" +  ent  + "</th>" 
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
	rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html)
	
	for lEntry in lDisplay_Output1:
		rc = Write_Rec_to_Html(4, '<tr>' , FileOut_Html)
		for ind in range(len(lEntry)):
			if lEntry[ind] == "Y":
				lEntry[ind] = " "
			sOut = "<td>" + str(lEntry[ind]) + "</td>"
			rc = Write_Rec_to_Html(8, sOut, FileOut_Html)
		rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '</table> ' , FileOut_Html)
	return (lDisplay_Output1) 
 
#******************************************************************************
#     Display EC Predictions                                                  *
#******************************************************************************
def Print_Predictions(AC,  dACs_Hits, FileOut_Html, lEC_Predictions, lNon_EC_Predictions): 
	lPredictions_Entry = [AC]
	
	tree = dACs_Hits[AC]["tree"]
	sHeader_Text = "Specific Peptide  Predictions:  " + dACs_Hits[AC] ["Description"]
	sHeader = '<h3 class="WhiteText">' + 	sHeader_Text + "</h3>"
	Write_Rec_to_Html(0,sHeader , FileOut_Html)
	
	rc = Write_Rec_to_Html(0, '<table border="1"> ' , FileOut_Html)
	
	rc = Write_Rec_to_Html(4, '<tr style="text-align: center;">' , FileOut_Html)
	lHeader = ["Function", "Prediction","Coverage"]
	for ent in lHeader:
		sOut = "<th>" +  ent  + "</th>" 
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
	rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html)
	
	for EC_Prediction in lEC_Predictions:
		rc = Write_Rec_to_Html(4, '<tr>' , FileOut_Html)
		rc = Write_Rec_to_Html(4, '<td>EC</td>' , FileOut_Html)
		sOut = "<td>" + EC_Prediction + "</td>" 
		EC_Prediction_Coverage = Calculate_EC_Prediction_Coverage(EC_Prediction, tree)
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
		sOut = "<td>" + str(EC_Prediction_Coverage)  + "</td>"
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
		rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html)
		
		lPredictions_Entry.append(["EC", EC_Prediction, EC_Prediction_Coverage])
	
	for Entry in  lNon_EC_Predictions:
		rc = Write_Rec_to_Html(4, '<tr>' , FileOut_Html)
		sOut = "<td>" + Entry[0] + "</td>"
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
		sOut = "<td>" + str(Entry[1]) + "</td>" 
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
		sOut = "<td>" + str(Entry[2]) + "</td>" 
		rc = Write_Rec_to_Html(4, sOut , FileOut_Html)
		rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html) 	
	
		lPredictions_Entry.append(Entry)
	
	rc = Write_Rec_to_Html(4, '</table> ' , FileOut_Html)
	
	lPredictions.append(lPredictions_Entry)
	return(0)

#******************************************************************************
#     Calculate The coverage for the predicted EC                             *
#****************************************************************************** 
def  Calculate_EC_Prediction_Coverage(EC_Prediction, tree):
	EC_Prediction_Coverage = 0
	sub_tree = tree.subtree(EC_Prediction)
	lSubTree_nids = [sub_tree[node].tag for node in  sub_tree.expand_tree(mode=Tree.DEPTH)]
	for nid in lSubTree_nids:
		node = sub_tree.get_node(nid)
		node_Coverage = node.data
		EC_Prediction_Coverage+=node_Coverage
	return(EC_Prediction_Coverage)
 
#******************************************************************************
#     Calculate SP Hits Distributions   (Empty plaxces don't count            *
#******************************************************************************
def Coverage(AC, dACs_Hits):
	dCoverage = dict()
	lProperties_Flat = reduce(operator.concat, dACs_Hits[AC]["Detailed_Properties_Map"])
	for iProperty in range(len(dACs_Hits[AC]["lProperties"])):
		dCoverage[iProperty] = lProperties_Flat.count(iProperty) 
	return  (dCoverage)
 
#******************************************************************************
#     Add  EC to the tree                                                     *
#******************************************************************************
def Add_EC_To_Tree(EC, AC, dACs_Hits,Coverage):
	EC = EC.replace(".-","")
	EC_Level = EC.count(".") + 1
 

	if  EC_Level == 1:
		ECLevel1 = EC.split(".")[0]
		if not tree.contains(ECLevel1):
			tree.create_node(ECLevel1, ECLevel1, parent = "Root" , data = Coverage)  
	if  EC_Level == 2:
		ECLevel1 = EC.split(".")[0]
		ECLevel2 = EC.split(".")[0]    + "." + EC.split(".")[1]  
		if not tree.contains(ECLevel1):
				tree.create_node(ECLevel1, ECLevel1, parent = "Root",  data = 0) 
		tree.create_node(ECLevel2, ECLevel2, parent = ECLevel1, data = Coverage ) 	

			
	if  EC_Level == 3:
		ECLevel1 = EC.split(".")[0]
		ECLevel2 = EC.split(".")[0]    + "." + EC.split(".")[1]      
		ECLevel3 = EC.split(".")[0]    + "." + EC.split(".")[1] + "." + EC.split(".")[2] 
		if not tree.contains(ECLevel1):
			tree.create_node(ECLevel1, ECLevel1, parent = "Root",    data = 0)
		if not tree.contains(ECLevel2):
			tree.create_node(ECLevel2, ECLevel2, parent = ECLevel1,    data = 0) 
		tree.create_node(ECLevel3, ECLevel3, parent = ECLevel2, data = Coverage )		
	
	if  EC_Level == 4:
		ECLevel1 = EC.split(".")[0]
		ECLevel2 = EC.split(".")[0]    + "." + EC.split(".")[1]      
		ECLevel3 = EC.split(".")[0]    + "." + EC.split(".")[1] + "." + EC.split(".")[2] 
		ECLevel4 = EC.split(".")[0]    + "." + EC.split(".")[1] + "." + EC.split(".")[2]  +  "." + EC.split(".")[3]
		if not tree.contains(ECLevel1):
			tree.create_node(ECLevel1, ECLevel1, parent = "Root",   data = 0)
		if not tree.contains(ECLevel2):
			tree.create_node(ECLevel2, ECLevel2, parent = ECLevel1, data = 0) 
		if not tree.contains(ECLevel3):
			tree.create_node(ECLevel3, ECLevel3, parent=ECLevel2,   data = 0)
		if not tree.contains(ECLevel4):
			tree.create_node(ECLevel4, ECLevel4, parent = ECLevel3, data = Coverage) 
	return(0)
	

 



 
#******************************************************************************
#     Predict EC from Tree                                                    *
#******************************************************************************
def Predict_EC_From_Tree(AC, dACs_Hits):
 
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
 
	 
 
#******************************************************************************
#     Break the sequence in chunks                                            *
#******************************************************************************
def chunkstring(string, length):
    return  (string[0+i:length+i] for i in range(0, len(string), length)) 

#******************************************************************************
#     Write HTML Records                                                      *
#******************************************************************************
def Write_Rec_to_Html(iPad, Rec, FileOut_Html):
		if fPrint_Html == False:
			return(0)
		sPad = " "*iPad
		sOutRec = sPad + Rec + cN
		FileOut_Html.write(sOutRec)
		return(0)

#******************************************************************************
#     Write the html file for AC                                              *
#******************************************************************************
def  Print_Html_For_AC (AC, lPrintOutput, dACs_Hits,  FileOut_Html, lHeader):
 
	if  len(dACs_Hits[AC]["Hits"] ) == 0:
		sHeader_Text = "Analysis of " + dACs_Hits[AC] ["Description"] + "  *** No hits found "
		sHeader = '<h5 class="White_Text">' + 	sHeader_Text + "</h5>"
		Write_Rec_to_Html(0,sHeader , FileOut_Html)
		return(0)
	
	sHeader_Text = "Analysis of " + dACs_Hits[AC] ["Description"]
	sHeader = '<h3 class="WhiteText">' + 	sHeader_Text + "</h3>"
	Write_Rec_to_Html(0,sHeader , FileOut_Html)
	rc = Write_Rec_to_Html(0, '<table border="1">' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '<tr style="text-align: center;">' , FileOut_Html)
	for ent in lHeader:
		sOut = "<th>" +  ent  + "</th>" 
		rc = Write_Rec_to_Html(8, sOut , FileOut_Html)
	rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html)
	
	for lOutput_Entry in lPrintOutput:
		if lOutput_Entry[2] != AC:
			continue
		rc = Write_Rec_to_Html(8, '<tr>' , FileOut_Html)
		for indx in range(len(lHeader)):
				my_attributes = ""
				sOutput = " "
				if  indx <  len(lOutput_Entry):
					sOutput =  lOutput_Entry[indx]
					
					if   lHeader[indx].startswith("Function")  and "Description"  not in lHeader[indx]:
						iColorNum = lSP_Types.index(lOutput_Entry[indx])
						sMyColor = lColor_Names_Singles[iColorNum+1]
						my_attributes = " class='" + sMyColor + "'"
				 
				sOut = "<td"  + my_attributes + ">" +  str(sOutput) + "</td>" 
				rc = Write_Rec_to_Html(12,sOut , FileOut_Html)
		rc = Write_Rec_to_Html(8,"</tr>" , FileOut_Html)
	
	rc = Write_Rec_to_Html(0,"</table>" , FileOut_Html)
	return(0)


 


#******************************************************************************
#     Build css                                                               *
#****************************************************************************** 
def  Build_Css(FileOut_Html, iTotal_Single_Colors):

	lBgColors = [ "White", "Gold", "Aqua","Pink", "Aquamarine",\
		"GreenYellow","Cornsilk",\
		"GoldenRod","Coral", "DarkMagenta", "Lawngreen", \
		"LawnGreen", "LightSkyBlue","PaleGOldenRod",\
		"Azure", "AntiqueWhite", "Cyan",\
		"Gainsboro", "Ivory",\
		"Lavender", "LightBlue", "LightSeaGreen"]
		
		
	lfontColors = ["Black","Blue", "Red"]
	lColor_Names_Singles = list()
	lColor_Names_Combinations = list()
	rc = Write_Rec_to_Html(0,"<style>" , FileOut_Html)
	ind = -1
	for bgColor in lBgColors:
		ind+=1
		sOut =  bgColor + "_Text"
		if ind < iTotal_Single_Colors: 
			lColor_Names_Singles.append(sOut)
		else:
			lColor_Names_Combinations.append(sOut)
		sOut = "." + sOut
		rc = Write_Rec_to_Html(4,sOut , FileOut_Html)
		rc = Write_Rec_to_Html(4,"{" , FileOut_Html)
		sOut = "color:" + "Black" + ";"
		rc = Write_Rec_to_Html(4,sOut , FileOut_Html)
		rc = Write_Rec_to_Html(4,"font-family: 'Courier New', monospace;" , FileOut_Html)
		sOut = 		"background-color: " + bgColor + ";"
		rc = Write_Rec_to_Html(4,sOut , FileOut_Html)
		rc = Write_Rec_to_Html(4,"}" , FileOut_Html)
	rc = Write_Rec_to_Html(0,"</style>" , FileOut_Html)
	return(lColor_Names_Singles, lColor_Names_Combinations )

#******************************************************************************
#     Determine      color                                                    *
#****************************************************************************** 
def   Determine_Color( Current_Attr, lColor_Names_Singles, lColor_Names_Combinations, lSeq_Pos_Attrs ):
		Current_Color = lColor_Names_Singles[0]
		if Current_Attr == []:
			Current_Color = lColor_Names_Singles[0]
		elif len(Current_Attr)  == 1:
			Current_Color = lColor_Names_Singles[Current_Attr[0]]
		elif len(Current_Attr)  == 2:
			lColorIndex = lSeq_Pos_Attrs[0]
			if  Current_Attr in lCombination_Colors_Used:
				iSelected_Color_Num = lCombination_Colors_Used.index(Current_Attr)
				Current_Color  = lColor_Names_Combinations[iSelected_Color_Num]
		return(Current_Color)

#******************************************************************************
#     Display Sequences                                                       *
#****************************************************************************** 
def Display_Sequence(AC, dACs_Hits ):
	if len(dACs_Hits[AC]["Hits"] ) == 0:
		return(0)

	
	sHeader = '<h3 class="WhiteText">' + 	dACs_Hits[AC] ["Description"] + "</h3>"
	Write_Rec_to_Html(0,sHeader , FileOut_Html)
	

	ind = -1
	sOut = "<p>"
	for sChunk in dACs_Hits[AC]["Chunks"]:
		ind+=1
		lSeq_Pos_Attrs = dACs_Hits[AC]["Seq_Positions_Attr"][ind]

		
		Prev_Attr = ["?"]
		for iPos in range(len(lSeq_Pos_Attrs)):
			Current_Attr = lSeq_Pos_Attrs[iPos]
			Next_Color = Determine_Color( Current_Attr, lColor_Names_Singles, lColor_Names_Combinations,lSeq_Pos_Attrs )
			if  Current_Attr != Prev_Attr:
				sOut = sOut +   '<span class="' + Next_Color + '">'
			sOut = sOut + sChunk[iPos]
			Prev_Attr = Current_Attr
		sOut = sOut + '</span>' + "<br>"  
	
	sOut = sOut + "</p>"
	rc = Write_Rec_to_Html(0,sOut , FileOut_Html)


	dACs_Hits[AC]["Coverage"] =  Coverage(AC,  dACs_Hits)
	return (0)
	
	
#******************************************************************************
#     Predict Function                                                        *
#****************************************************************************** 
def Load_EC_Tree(AC, dACs_Hits):
	dACs_Hits[AC]["Predictions"] = list()
	for ent in dACs_Hits[AC]["Functions_Coverage"]:
		if ent[0] != "EC":
			dACs_Hits[AC]["Predictions"].append(ent)
			continue
		Coverage = ent[2]
		if Coverage < dACs_Hits[AC]["SP_Len_Thresh"]:
			continue
		EC = ent[1]
		rc = Add_EC_To_Tree(EC, AC, dACs_Hits,Coverage)
	return(0)
  
#******************************************************************************
#     Print Legend                                                            *
#****************************************************************************** 
def   Print_Legend( AC,  dACs_Hits, FileOut_Html,  lColor_Names_Singles, lColor_Names_Combinations,  lSP_Types   ):
	lSingle_Attrs = list()
	lComb_Attrs = list()
	for lAttrs in dACs_Hits[AC]["Seq_Positions_Attr"] :
		for sAttr in lAttrs:
			if len(sAttr) == 0:
				continue

			if len(sAttr ) == 1 :
				lSingle_Attrs.append(sAttr[0])
			if len(sAttr ) > 1:
				lComb_Attrs.append(sAttr)			

	lSingle_Attrs = list(set(lSingle_Attrs))
		

	if len(lComb_Attrs) > 0:
		lComb_Attrs_New = list()
		dedup = map(list, unique(map(tuple, lComb_Attrs)))
		for ent in dedup:
			lComb_Attrs_New.append(ent)
		lComb_Attrs = lComb_Attrs_New
	
	rc = Write_Rec_to_Html(0, '<h4>Legend</h4>' , FileOut_Html)
	rc = Write_Rec_to_Html(0, '<table border="1">' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '<tr style="text-align: center;">' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '<th>SP Type</th>' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '<th>Color</th>' , FileOut_Html)
	rc = Write_Rec_to_Html(4, '</tr>' , FileOut_Html)
	
	for ind_Color in lSingle_Attrs:
		if ind_Color == 0:
			continue
		rc = Write_Rec_to_Html(8, '<tr>' , FileOut_Html)
		SP_Type = lSP_Types[ind_Color -1]
		sOut = '<td>' +  SP_Type + '</td>' 
		rc = Write_Rec_to_Html(12, sOut , FileOut_Html)

		sColorName = lColor_Names_Singles[ind_Color]
		sOut = "<td class='" + sColorName + "'>" + "  " +  "</td>" 
		rc = Write_Rec_to_Html(12, sOut , FileOut_Html)
		rc = Write_Rec_to_Html(8, '</tr>' , FileOut_Html)
		
	for eComb in lComb_Attrs:
		iSP_Type1_Num = eComb[0]
		iSP_Type2_Num = eComb[1]
		rc = Write_Rec_to_Html(8, '<tr>' , FileOut_Html)
		sOut = '<td>' + " Overlap between  " +  lSP_Types[iSP_Type1_Num - 1] + " and " + lSP_Types[iSP_Type2_Num -1] + '</td>' 
		rc = Write_Rec_to_Html(12, sOut , FileOut_Html)
		if eComb in lComb_Attrs:
			iColor = lComb_Attrs.index(eComb)
			sColorName1 = lColor_Names_Combinations[iColor]
		else:
			sColorName1 = lColor_Names_Singles[0]
		sOut = "<td class='" + sColorName1 + "'>" + " " +  "</td>" 
		rc = Write_Rec_to_Html(12, sOut , FileOut_Html)
		rc = Write_Rec_to_Html(8, '</tr>' , FileOut_Html)
		
	
	rc = Write_Rec_to_Html(4, '</table> ' , FileOut_Html)
 
	return (0)
	

	
	

#******************************************************************************
#     Program Starts here                                                     *
#****************************************************************************** 
StartTime = datetime.datetime.now()

CommonArea = read_params( sys.argv )   
parser = CommonArea['parser'] 
results = parser.parse_args()
dParms  = Prepare_Run_Parameters(results)

 



print("Creating html analysis html output into ", dParms["Output_Html_FileName"])
if   dParms["Output_Pred_xlsx_Required"].lower().startswith("y"):
	print("Creating predictions excel summary into  ", dParms["Output_Predictions_xlsx_FileName"])

fPrint_Html = True
if dParms["Print_Html_Report"].lower().startswith("n"):
	fPrint_Html = False
	FileOut_Html = None
	print("** Html Analysis report will not be printed **")
	
if fPrint_Html == True:
	Output_Html_FileName = dParms["Output_Html_FileName"]
	FileOut_Html =  open(Output_Html_FileName,"w")

iDisplayAfter = 100
cntSeq = 0
iLenChunk = dParms["w_Seq_Display"]
 
dHitMap = dict() 
lOutput = list()
lPrintOutput = list()
lChunks = list()
dACs_Hits = dict()
cN ="\n"
lAC_Acttributes = list()
lHeader = ["Hit_Location_Start", "Hit_Location_End", "Protein_ID", "SP", \
			"Function",  "Function_Description"]
iTotal_Single_Colors = 7
lSingle_Colors_Used = list()
lCombination_Colors_Used = list()
lSP_Types = ["EC","ZF", "GPCR_OR", "GPCR_NOR"]
lProcessed_ACs = list()
lPredictions = list() 


lColor_Names_Singles, lColor_Names_Combinations  = Build_Css(FileOut_Html, iTotal_Single_Colors)

with open(dParms["dSPs_File_Name"]) as json_file:
	dSPs = json.load(json_file)
Input_Fasta_File =   dParms["Input_Fastq_Name"]
 
	
myfasta_sequences = SeqIO.parse(open(Input_Fasta_File),'fasta')
for seq in myfasta_sequences:
	Seq = str(seq.seq)
	iLenSeq = len(Seq)
	AC = seq.id
	if AC in lProcessed_ACs:
		print("AC: ", AC, " has been processed,  found duplicate in input and ignored")
		continue
	lProcessed_ACs.append(AC)
	dACs_Hits[AC] = dict()
	dACs_Hits[AC]["Description"] = seq.description
	dACs_Hits[AC]["Hits"] = list()
	dACs_Hits[AC]["Attributes"] = list()
	lChunks =   list(chunkstring(Seq, iLenChunk) )
	dACs_Hits[AC]["Chunks"] = lChunks
	dACs_Hits[AC]["Seq"] = Seq
	dACs_Hits[AC]["SeqLen"] = len(Seq)
	dACs_Hits[AC]["Seq_Attr"] =  [ [] for _ in range(dACs_Hits[AC]["SeqLen"]) ]
	dACs_Hits[AC]["SP_Len_Thresh"] = dParms["SP_Len_Thresh"]
	SP_Len_Thresh = dACs_Hits[AC]["SP_Len_Thresh"]

	
	cntSeq+=1
	if cntSeq % iDisplayAfter == 0:
		Current_Time  = datetime.datetime.now() 
		print("Read and searched SPs in ", str(cntSeq), "  Sequences ",  Current_Time.strftime("%Y-%m-%d %H:%M:%S"))
 
	for Motif in dSPs.keys():
		##if len(Motif) < SP_Len_Thresh:
			###continue
		iCntHits = Seq.count(Motif)

		if iCntHits > 0:
			lHits_Locations0 = [m.start() for m in re.finditer(Motif, Seq)]
			lHits_Locations = [e+1 for e in lHits_Locations0]
			lHitEntry = [Motif,  lHits_Locations]
			dACs_Hits[AC]["Hits"].append(lHitEntry)
	
for AC in sorted(dACs_Hits.keys()):
	tree = None
	tree = Tree()
	tree.create_node("Root", "Root" , data = 0)  # root node
	
	if len(dACs_Hits[AC]["Hits"]) == 0:
		continue
	dACs_Hits[AC]["Properties_Map"] = [ [] for _ in range(dACs_Hits[AC]["SeqLen"]) ]
	seqLen = dACs_Hits[AC]["SeqLen"]
	dACs_Hits[AC]["Detailed_Properties_Map"] = [ [] for _ in range(dACs_Hits[AC]["SeqLen"]) ]
	dACs_Hits[AC]["lProperties"] = list()
	SeqLen = len(dACs_Hits [AC]["Seq"])
	lOutput  =  list()
	for lHit_Entry in dACs_Hits[AC]["Hits"]:
		SP = lHit_Entry[0]
		SP_Property = dSPs[SP]
		lHits_Location = lHit_Entry[1]
		iCntHits_SP_On_Protein = len(lHits_Location)
		for ind in range(iCntHits_SP_On_Protein):
			lMap_Hit_Entry = [lHits_Location[ind],lHits_Location[ind] + len(SP) -1, AC,SP]

			for Fn in sorted(dSPs[SP].keys()):
				Fn_Description = dSPs[SP][Fn]
					
				if Fn == "ZF" and Fn_Description == "Y":
					Fn_Description = " "
				if Fn == "GPCR_OR" and Fn_Description == "Y":
					Fn_Description = " "
				lMap_Hit_Entry.append(Fn)
				lMap_Hit_Entry.append(Fn_Description)
				

####    =============>    lSP_Types = ["EC","ZF", "GPCR_OR", "GPCR_NOR"]
				for indz in range(lHits_Location[ind] -1 ,lHits_Location[ind] + len(SP) -1):
						for indy in range(len(lSP_Types)):
							if Fn == lSP_Types[indy]:
								dACs_Hits[AC]["Seq_Attr"][indz].append(indy + 1)
								dProperty = dSPs[SP]
								if dProperty not in dACs_Hits[AC]["lProperties"]:
									dACs_Hits[AC]["lProperties"].append(dProperty)
								iProperty = dACs_Hits[AC]["lProperties"].index(dProperty)
								dACs_Hits[AC]["Detailed_Properties_Map"][indz].append(iProperty)
								break

			lOutput.append(lMap_Hit_Entry)
	lOutput = sorted(lOutput)
	for indx in range(len(lOutput)):
		lPrintOutput.append(lOutput[indx])
	
	dACs_Hits[AC]["Seq_Positions_Attr0"] = list()
	for lSeq_Pos_Attr in dACs_Hits[AC]["Seq_Attr"]:
		lSeq_Pos_Attr1 =  sorted(list(set(lSeq_Pos_Attr)))
		if len(lSeq_Pos_Attr1) == 0:
			if  0 not in lSingle_Colors_Used:
				lSingle_Colors_Used.append(0)
		elif len(lSeq_Pos_Attr1) == 1:
			if  lSeq_Pos_Attr1[0] not in lSingle_Colors_Used:
				lSingle_Colors_Used.append(lSeq_Pos_Attr1[0])
		elif len(lSeq_Pos_Attr1) == 2:
			if lSeq_Pos_Attr1 not in lCombination_Colors_Used:
				lCombination_Colors_Used.append(lSeq_Pos_Attr1)
	
		dACs_Hits[AC]["Seq_Positions_Attr0"].append( lSeq_Pos_Attr1)
	n = iLenChunk
	dACs_Hits[AC]["Seq_Positions_Attr"] = [dACs_Hits[AC]["Seq_Positions_Attr0"][i:i+n] for i in range(0, len(dACs_Hits[AC]["Seq_Positions_Attr0"]), n)]
	dACs_Hits[AC]["Seq_Positions_Attr_Map"] = dACs_Hits[AC]["Seq_Positions_Attr0"]
	dACs_Hits[AC]["Seq_Attr"]
	del dACs_Hits[AC]["Seq_Positions_Attr0"]
	del dACs_Hits[AC]["Seq_Attr"]
	
	
for AC in sorted(dACs_Hits.keys()):
	rc = Print_Html_For_AC (AC, lPrintOutput, dACs_Hits,  FileOut_Html, lHeader) 
	rc = Display_Sequence(AC, dACs_Hits )
	if len(dACs_Hits[AC]["Hits"]) > 0:
		rc = Print_Legend(AC,  dACs_Hits,  FileOut_Html,   lColor_Names_Singles, lColor_Names_Combinations, lSP_Types   )
		dACs_Hits[AC]["Functions_Coverage"]  = Display_Coverage(AC, dACs_Hits,  FileOut_Html) 
		tree = None
		tree = Tree()
		tree.create_node("Root", "Root" , data = 0)  # root node

		rc = Load_EC_Tree(AC, dACs_Hits )
		dACs_Hits[AC]["tree"] = tree
		lEC_Predictions  = Predict_EC_From_Tree(AC, dACs_Hits)

		lNon_EC_Predictions = list()

		for lEntry in dACs_Hits[AC]["Functions_Coverage"]:
			if lEntry[0] == "EC":
				continue

			lDisplay_Function = lEntry
			lNon_EC_Predictions.append(lDisplay_Function)
			
		if len(lEC_Predictions) > 0  or len(lNon_EC_Predictions)  > 0 :
			rc = Print_Predictions(AC,  dACs_Hits, FileOut_Html, lEC_Predictions, lNon_EC_Predictions) 


	del dACs_Hits[AC]
	
lPrintOutput.insert(0, lHeader)

lSingle_Colors_Used = sorted(lSingle_Colors_Used)
lCombination_Colors_Used = sorted(lCombination_Colors_Used)



if dParms["Output_xlsx_Required"].lower().startswith("y") \
or  dParms["Output_Simple_Html_Required"].lower().startswith("y"):
		df = pd.DataFrame(lPrintOutput)

if dParms["Output_xlsx_Required"].lower().startswith("y"): 
	print("Creating excel output into ", dParms["Output_xlsx_FileName"])
	df.to_excel(dParms["Output_xlsx_FileName"])

if dParms["Output_Simple_Html_Required"].lower().startswith("y"):
	print("Creating simple html output into ", dParms["Output_Simple_Html_FileName"])
	df.to_html(dParms["Output_Simple_Html_FileName"]) 

 
if   dParms["Output_Pred_xlsx_Required"].lower().startswith("y"):
	df_preds = pd.DataFrame(lPredictions)
	df_preds.to_excel(dParms["Output_Predictions_xlsx_FileName"])


if fPrint_Html == True:
	FileOut_Html.close()
	
EndTime = datetime.datetime.now()
print ("Start Time: ", StartTime.strftime("%Y-%m-%d %H:%M:%S"))
print ("End Time: ", EndTime.strftime("%Y-%m-%d %H:%M:%S"))
print("Program Completed Successfully")
sys.exit(0)