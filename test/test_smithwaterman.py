from hw3align.sequences import *
from hw3align.getmatrix import blosum50
from hw3align.smithwaterman import sw

#As a first test, lets make sure that the trace and alignment work okay regardless of starting 
#or ending with a gap: 
def sw_Gap_beginB_test():
	assert sw("ACDEFG", "EFG", blosum50, -3, -1)[3:] == ("EFG", "EFG") 
	#recall that sw() takes two strings to align, a scoring matrix, and gap start/extend penalties
def sw_Gap_beginA_test():
	assert smithWaterman("EFG", "ACDEFG", blosum50, -3, -1)[3:] == ("EFG", "EFG")
	#this is like above but now seq a has the gap
#def test_SW_endAGap():
	#assert smithWaterman("ACD", "ACDEFG", blosum50, -3, -1)[3:] == ("ACD", "ACD")
#def test_SW_endBGap():
	#assert smithWaterman("ACDEFG", "ACD", blosum50, -3, -1)[3:] == ("ACD", "ACD")

# I want to make sure I can trace through matches and mismatches:
def test_SW_trace_through_mismatch():
	return None 
	#assert smithWaterman("ACDAFG", "ACDEFG", blosum50, -3, -1)[2:] == (41.0, 'ACDAFG', 'ACDEFG')
# I want to make sure I can trace through gaps and extensions in both strings:
#def test_SW_trace_through_indel_B():
	#assert smithWaterman("GARRETT", "GAETT", blosum50, -3, -1)[3:] == ('GARRETT', 'GA--ETT')
#def test_SW_trace_through_indel_A():
	#assert smithWaterman("GAETT", "GARRETT", blosum50, -3, -1)[3:] == ('GA--ETT', 'GARRETT')

# I pulled a random pair and ran it through EMBOSS
# to get the optimal alignment to check.
# If it sticks around a while, the results should be at
# https://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=emboss_water-I20190223-011138-0845-38663360-p1m
#testA, testB = posPairs[5]
# This pair are the sequences 
#('RFKWGPASQQILFQAYERQKNPSKEERETLVEECNRAECIQRGVSPSQAQGLGSNLVTEVRVYNWFANRRKEEAFRH', 'GRKRKIDRDAVLNMWQQGLGASHISKTMNIARSTVYKVINESN')
#def test_SW_realData():
	#start, end, bestScore, alignA, alignB = smithWaterman(testA, testB, blosum50, -10, -1)
	#assert start == (35, 1)
	#assert end == (59, 26)
	#assert bestScore == 35
	#assert alignA == "RAECIQR-GVSPSQAQGLGSNLVTE"
	#assert alignB == "RKRKIDRDAVLNMWQQGLGASHISK"
