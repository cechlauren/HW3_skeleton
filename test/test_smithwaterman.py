from hw3align.sequences import *
from hw3align.getmatrix import blosum50
from hw3align.smithwaterman import sw
#############################################################################################
#As a first test, lets make sure that the trace and alignment work okay regardless of starting 
#with a gap: 
def test_sw_Gap_beginB():
	assert sw("EFGACD", "ACD", blosum50, -3, -1)[3:] == ("ACD", "ACD") 
	#recall that sw() takes two strings to align, a scoring matrix, and gap start/extend penalties
def test_sw_Gap_beginA():
	assert sw("ACD", "EFGACD", blosum50, -3, -1)[3:] == ("ACD", "ACD")
	#this is like above but now seq a has the gap
#############################################################################################
#As a second test, lets make sure that the trace and alignment work okay regardless of ending
#with a gap: 
def test_sw_Gap_endA():
	#assert sw("ACD", "ACDEFG", blosum50, -3, -1)[3:] == ("ACD", "ACD")
def test_sw_Gap_endB():
	#assert sw("ACDEFG", "ACD", blosum50, -3, -1)[3:] == ("ACD", "ACD")
#############################################################################################
#As a third test, lets make sure sw() can trace thru matches & mismatches:
def test_sw_trace_mismatch():
	assert sw("ACDAFG", "ACDEFG", blosum50, -3, -1)[2:] == (41.0, 'ACDAFG', 'ACDEFG') #make sure the score is right
#############################################################################################	
#As a final test, make sure sw() can trace through gaps and extensions in both strings:
def test_sw_trace_indelB():
	assert sw("LAREN", "LARN", blosum50, -3, -1)[3:] == ('LAREN', 'LAR-N')
def test_sw_trace_indelA():
	assert sw("LARN", "LAREN", blosum50, -3, -1)[3:] == ('LAR-N', 'LAREN')


