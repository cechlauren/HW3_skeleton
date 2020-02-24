import sys
from hw3align.sequences import *
from hw3align.optimalmatrix import *
from hw3align.getmatrix import blosum50
#######################################################################################################
#To test the optimization of these matrices, need to get true positive and true negative alignments

TP_Alignments, TN_Alignments = getAlignments(scorematrix=blosum50, gap_start=-10, gap_extend= -1)

def test_get_alignments():
	assert TP_Alignments[7] == ("NSNQIKILGNQGSFLTKG-PSKLNDRADSRRSLW--------DQGNFPLIIK------NLKI", "NCSTFYVVKEDGTIVYTGTATSMFD-NDTKETVYIADFSSVNEEGTYYLAVPGVGKSVNFKI")
#making sure the alignment example in my sw testing looks as expected.   

#######################################################################################################
#Now the scoring from sw will be redone

def test_scoring():
	
    testA, testB = TP_Alignments[7]
    assert scoreAlignment(testA, testB, blosum50, -10, -1) == 35.0
#making sure that the scoring example in my sw testing looks as expected
#######################################################################################################
#Utilize running output to verify that actual function is running correctly.
#Confirm data stays the right length:

def test_geneticAlg():
	#return None
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		blosum50, 1, 1, 1, 10, 5, 5, 3, -10, -1, true_pos_align, true_neg_align)
	assert len(pop) == 10
	#assert len(scores) == 10
	#assert len(library) == 3
	#assert len(objectiveMeans) == 5 + 1 # there's an extra -inf at the beginning so that it starts off increasing

#######################################################################################################
#If time, will test objective function:
#true_pos_aligns, true_neg_aligns = getAlignments(scorematrix=blosum50, gap_start=-10, gap_extend=-1)
#def test_objective_function():
	#recall Obj. fxn: def scoreMatrix(true_pos_align, true_neg_align, scorematrix, gap_start, gap_extend):
	#true_pos_align = true_pos_aligns[7, ?, ?]
	#true_neg_align = true_neg_aligns[?, ?, ?]
	#assert scoreMatrix(true_pos_align, true_neg_align, blosum50, -10, -1) == XXX

#######################################################################################################	
#Also important to check that the matrix changes but will not make excessive # of copies
#def test_mutate_matrix():
	#return None
    #M = blosum.50.copy()
    #M2 = mutateMatrix(M, 1, 1)
    #assert M2 = M
    #assert M2 != blosum50
    

#def test_selection():
	#return None
	# Very random and hard to test, but I want to make sure that my weights
	# are working and that I'm sampling with replacement: 
	#L = ["a", "b", "c"]
	#w = [0, 0, 1]
	#assert selection(L, w) == ["c", "c", "c"]

#def test_scale_scores():
	#return None
	
	# just make sure I'm doing what I said: largest/smallest is 10^selectivePressure
	# and the sum is 1.
	#a,b = scaleScores([1, 2], 1)
	#assert a+b == 1
	#assert b/a == 10
	


		



















