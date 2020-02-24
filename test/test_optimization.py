import sys
from BMI203_HW3_alignment.sequences import *
from BMI203_HW3_alignment.optimizeMatrix import *
from BMI203_HW3_alignment.matrices import *
import numpy as np
from smith_waterman import algs

tpAlignments, tnAlignments = getAlignments(scoringMatrix=blosum50, gapStart=-10, gapExtend= -1)

def test_get_alignments():
	return None
    #assert tpAlignments[5] == ("RAECIQR-GVSPSQAQGLGSNLVTE", "RKRKIDRDAVLNMWQQGLGASHISK")
   


#will redo scoring from sw
def test_scoring():
	return None
    #testA, testB = tpAlignments[5]
    #assert scoreAlignment(testA, testB, blosum50, -10, -1) ==35
    
#no test of obj fxn code since lots of OH.

#make sure matrix changes but wont make more copies than necessary
def test_mutate_matrix():
	return None
    #M = blosum.50.copy()
    #M2 = mutateMatrix(M, 1, 1)
    #assert M2 = M
    #assert M2 != blosum50
    

def test_selection():
	return None
	# Very random and hard to test, but I want to make sure that my weights
	# are working and that I'm sampling with replacement: 
	#L = ["a", "b", "c"]
	#w = [0, 0, 1]
	#assert selection(L, w) == ["c", "c", "c"]

def test_scale_scores():
	return None
	
	# just make sure I'm doing what I said: largest/smallest is 10^selectivePressure
	# and the sum is 1.
	#a,b = scaleScores([1, 2], 1)
	#assert a+b == 1
	#assert b/a == 10
	

# For the actual meat of the function, I mostly can only use the running output to
# verify that it's running correctly. Here, I make sure that all our data stays
# the right length:
def test_geneticAlg_skeleton():
	return None
	#pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		#blosum50, 1, 1, 1, 10, 5, 5, 3, -10, -1, tpAlignments, tnAlignments)
	#assert len(pop) == 10
	#assert len(scores) == 10
	#assert len(library) == 3
	#assert len(objectiveMeans) == 5 + 1 # there's an extra -inf at the beginning so that it
										# starts off increasing
		



















