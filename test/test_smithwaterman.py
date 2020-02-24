from hw3align.sequences import *
from hw3align.getmatrix import blosum50
from hw3align.smithwaterman import sw
###############################################################################################
#As a first test, lets make sure that the trace and alignment work okay regardless of starting 
#with a gap: 
def sw_Gap_beginB_test():
	assert sw("EFGACD", "ACD", blosum50, -3, -1)[3:] == ("ACD", "ACD") 
	#recall that sw() takes two strings to align, a scoring matrix, and gap start/extend penalties
def sw_Gap_beginA_test():
	assert sw("ACD", "EFGACD", blosum50, -3, -1)[3:] == ("ACD", "ACD")
	#this is like above but now seq a has the gap
#############################################################################################
#As a second test, lets make sure that the trace and alignment work okay regardless of ending
#with a gap: 
def sw_Gap_endA_test():
	assert sw("ACD", "ACDEFG", blosum50, -3, -1)[3:] == ("ACD", "ACD")
def sw_Gap_endB_test():
	assert sw("ACDEFG", "ACD", blosum50, -3, -1)[3:] == ("ACD", "ACD")
#############################################################################################
#As a third test, lets make sure sw() can trace thru matches & mismatches:
def sw_trace_mismatch_test():
	assert sw("ACDAFG", "ACDEFG", blosum50, -3, -1)[2:] == (41.0, 'ACDAFG', 'ACDEFG') #make sure the score is right
#############################################################################################	
#As a final test, make sure sw() can trace through gaps and extensions in both strings:
def sw_trace_indelB_test():
	assert sw("LARENCECH", "LARNCECH", blosum50, -3, -1)[3:] == ('LARENCECH', 'LAR-NCECH')
def sw_trace_indelA_test():
	assert sw("LARNCECH", "LARENCECH", blosum50, -3, -1)[3:] == ('LAR-NCECH', 'LARENCECH')


