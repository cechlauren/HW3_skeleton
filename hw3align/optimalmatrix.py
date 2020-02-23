import sys
import random
import numpy as np
from .sequences import *
from .smithwaterman import sw

def getAlignments(scorematrix, gap_start, gap_extend):
	'''
	This function will retrieve the true positive and true negative alignments using the scoring matrix, gap_start,
	and gap_extend costs.
	This function will return a list of true positive alignments (tuples of strings) and a list of true negative alignments
	'''
	true_pos_align = [] #list
	for a,b in positives: #the positive pairs we identified earlier from the positive pairs txt file
		score, matchedA, matchedB = sw(a,b, scorematrix, gap_start, gap_extend)[2:]
		true_pos_align.append((matchedA, matchedB)) #add those matched pairs to the list
	true_neg_align = []
	for a,b in negatives: #the negative pairs we identified earlier from the negative pairs txt file
		score, matchedA, matchedB = sw(a,b, scorematrix, gap_start, gap_extend)[2:]
		true_neg_align.append((matchedA, matchedB)) #add those negative pairs to the list
	return true_pos_align, true_neg_align

def scoreAlignment(a, b, scorematrix, gap_start, gap_extend):
	'''
	This function will take an aligmment, score it and return that score.
	The input for this function are aligned strings (a and b),
	a scoring matrix, and the gap start/extend parameters.
	
	'''
	openGapA = False #start off with no gaps
	openGapB = False
	score = 0 #start of with no score
	assert len(a) == len(b) #alignments need to be the same length
	for A, B in zip(a, b):
		if A == "-" and openGapA == False:
			score += gap_start #needed to add the gapstart cost
			openGapA = True
		elif A == "-" and openGapA == True: #if there was alread a gap start
			score += gap_extend #then we do extension if still no match
		elif B == "-" and openGapB == False:
			score += gap_start
			openGapB = True
		elif B == "-" and openGapB == True:
			score += gap_extend
		else:
			score += scorematrix[frozenset((A, B))] #for all other cases need to ref the score matrix
			openGapA = False
			openGapB = False
	return score

scoreMatrixD = {} #global dictionary
def scoreMatrix(true_pos_align, true_neg_align, scorematrix, gap_start, gap_extend):
	'''
	This function will score the positive or neg. alignments using the scorematrix.
	The objective fxn is the SUM of true positive rates at false positive rates of 0.0, 0.1, 0.2, and 0.3.
	The function's inputs are true positive and true negative TUPLES of aligned STRINGS, a scoring matrix,
	and the gap start/extend costs for alignment of those strings.
	The function returns the objective function.
	'''
	#Will first determine the cutoff for each false positive rate.
	k = frozenset(scorematrix.items()) #make the scoring matrix an immutable set of values that reference as k
	if k in scoreMatrixD.keys(): #if there is a match then
		return scoreMatrixD[k]
	negScores = [scoreAlignment(a,b, scorematrix, gap_start, gap_extend) for a,b in true_neg_align] #the bad scores
	cutoffs = np.percentile(negScores, [100, 90, 80, 70]) #designate some cutoffs
	#Then find the true positve rates for each cutoff
	posScores = [scoreAlignment(a,b, scorematrix, gap_start, gap_extend) for a,b in true_pos_align]
	#The best score will be the one where the average positive scores pass the cutoff based on the FPR we designated
	score = sum([sum(map(lambda x: x > cutoff, posScores))/len(posScores) for cutoff in cutoffs])
	scoreMatrixD[k] = score
	return score

def mutateMatrix(scorematrix, mutationChance, mutationAmount):
	'''
	This function will mutate values of the matrix where each value has a [0,1] chance of mutating,
	and will change by adding some gaussian noise with a mean of 0 and standard deviation that we call
	mutationAmount.
	mutationChance= probability each entry in a scoring matrix will mutate
	mutationAmount= standard deviation of mutation in gaussian dist.
	The function returns the same matrix, just modified.
	'''
	for key, value in scorematrix.items():
		if random.random() < mutationChance: #if some random value is less than the chance of mutation then
			scorematrix[key] = value + random.gauss(0, mutationAmount) #new key has that value and some random noise
	return scorematrix

def selection(pop, weights):
	'''
	This function will create a a whole new kind of matrices by sampling with replacement and with weights.
	The function input is a LIST of objects and a LIST of selection probabilities.
	The function returns a LIST of the same size.
	'''
	return list(np.random.choice(pop, size=len(pop), p=weights))

def scaleScores(scores, selectivePressure):
	'''
	This function will scale scores so that the largest is 10^selectivePressure and the  smallest sum is 1.
	Selective Pressure describes the most fit individual scaled to 10^selectivePressure * least fit
	The function input is a list of scores and returns a new list of scaled scores.
	'''
	oldMin = min(scores) #old minimum value of the scores
	oldMax = max(scores) #old maximum value of the scores
	newMin = 1 #new sum is 1
	newMax = 10**selectivePressure #selective pressure for that amino substitution
	newScores =  [(newMax - newMin)*(x - oldMin)/(oldMax - oldMin) + newMin for x in scores]
	#the new scores report the difference between new selective pressure times difference of scores and old min 
	#scores divided by difference in old scores
	return [s/sum(newScores) for s in newScores]

#To optimize a matrix, utilize a genetic algorithm
def optimizeMatrix_geneticAlg(scorematrix, mutationChance, mutationAmount,
	selectivePressure, N, totalStepsToStop, stepsWithNoImprovement,
	librarySize, gap_start, gap_extend, true_pos_align, true_neg_align):
	'''
	This function will optimize an alignment score matrix using a genetic algorithm.
	The function stops if it hits "totalItersToStop" iterations, or if it doesn't see
	a new objective function value in the top (librarySize) in stepsWithNoImprovement steps.
	The input of this function is the scorematrix, mutationChance, mutationAmount,		
	selectivePressure, n= population size, 
	"totalStepsToStop" =our runtime step cutoff, 
	"stepsWithNoImprovement"= the step cutoff for seeing no new entries in library of best matrices,
	librarySize= "number of best matrices to remember and re-seed population with",			
	gap_start, gap_extend, true_pos_align, and finally true_neg_align.
	
	The function will return the final population of matrices, the scores for matrices at the last step taken,
	the library of best matrices, and a list of the mean objective function value at each step taken.
	'''
	# initialize population
	pop = [scorematrix.copy() for i in range(N)]
	library = {}
	libraryMin = float("-inf")
	# loop, keeping track of change in objective function
	objectiveMeans = []
	noImprovementCounter = 0
	while True:
		noImprovementCounter += 1
		# mutate
		pop = [mutateMatrix(m, mutationChance, mutationAmount) for m in pop]
		# score
		scores = [scoreMatrix(true_pos_align, true_neg_align, m, gap_start, gap_extend) for m in pop]
		# bank good ones
		for score, mat in zip(scores, pop):
			if score > libraryMin:
				if score not in library.keys():
					noImprovementCounter = 0
					library[score] = mat.copy()
					if len(library) > librarySize:
						libraryMin = min(library.keys())
						del library[libraryMin]
						libraryMin = min(library.keys())
		# assess.
		objectiveMean = sum(scores)/N
		# print (len(objectiveMeans), objectiveMean, max(scores), max(library.keys()))
		objectiveMeans.append(objectiveMean)
		# Should we stop?
		if len(objectiveMeans) > totalStepsToStop or noImprovementCounter > stepsWithNoImprovement:
			break
		# we're continuing. Let's make a new generation.
		# Rescale objective to get weights
		weights = scaleScores(scores, selectivePressure)
		# select a new generation
		pop = selection(pop, weights)[:-len(library)] + [x.copy() for x in library.values()]
		# seed with good ones
	return pop, scores, library, objectiveMeans
