import sys
import random
import numpy as np
from .sequences import *
from .smithwaterman import sw

def getAlignments(scorematrix, gap_start, gap_extend):
	'''
	fetches true positive and true negative alignments using the scoring matrix, gap_start,
	and gap_extend costs.
	returns a list of true positive alignments (tuples of strings) and the same for true negatives.
	'''
	true_pos_align = []
	for a,b in positives:
		score, matchedA, matchedB = sw(a,b, scorematrix, gap_start, gap_extend)[2:]
		true_pos_align.append((matchedA, matchedB))
	true_neg_align = []
	for a,b in negatives:
		score, matchedA, matchedB = sw(a,b, scorematrix, gap_start, gap_extend)[2:]
		true_neg_align.append((matchedA, matchedB))
	return true_pos_align, true_neg_align

def scoreAlignment(a, b, scorematrix, gap_start, gap_extend):
	'''
	given an alignment, just score it.
	takes as input a and b; aligned strings,
	as well as a scoring matrix and gap/extend parameters.
	returns a score for this alignment.
	'''
	openGapA = False
	openGapB = False
	score = 0
	assert len(a) == len(b)
	for A, B in zip(a, b):
		if A == "-" and openGapA == False:
			score += gap_start
			openGapA = True
		elif A == "-" and openGapA == True:
			score += gap_extend
		elif B == "-" and openGapB == False:
			score += gap_start
			openGapB = True
		elif B == "-" and openGapB == True:
			score += gap_extend
		else:
			score += scorematrix[frozenset((A, B))]
			openGapA = False
			openGapB = False
	return score

scoreMatrixD = {}
def scoreMatrix(true_pos_align, true_neg_align, scorematrix, gap_start, gap_extend):
	'''
	score our alignments using scorematrix.
	our objective function is the sum of TP rates at FP rates of 0.0,
	0.1, 0.2, and 0.3.
	takes as input true positive and true negative tuples of aligned strings, as
	well as the scoring matrix and gap/extend costs for alignment.
	returns the objective function.
	'''
	# First, find the cutoff for each FP rate.
	k = frozenset(scorematrix.items())
	if k in scoreMatrixD.keys():
		return scoreMatrixD[k]
	negScores = [scoreAlignment(a,b, scorematrix, gap_start, gap_extend) for a,b in true_neg_align]
	cutoffs = np.percentile(negScores, [100, 90, 80, 70])
	# then, find the TP rates for each cutoff.
	posScores = [scoreAlignment(a,b, scorematrix, gap_start, gap_extend) for a,b in true_pos_align]
	score = sum([sum(map(lambda x: x > cutoff, posScores))/len(posScores) for cutoff in cutoffs])
	scoreMatrixD[k] = score
	return score

def mutateMatrix(scorematrix, mutationChance, mutationAmount):
	'''
	mutate some of the values of the matrix.
	each value has mutationChance probability [0, 1] of mutating,
	and will change by adding gaussian noise with mean 0 and standard deviation mutationAmount.
	returns the SAME matrix, modified.
	'''
	for key, value in scorematrix.items():
		if random.random() < mutationChance:
			scorematrix[key] = value + random.gauss(0, mutationAmount)
	return scorematrix

def selection(pop, weights):
	'''
	creates a new generation of matrices by sampling with replacement with weights.
	takes as input a list of objects and a list of selection probabilities, and
	returns a listof the same size.
	'''
	return list(np.random.choice(pop, size=len(pop), p=weights))

def scaleScores(scores, selectivePressure):
	'''
	scales scores so the largest/smallest is 10^selectivePressure and the sum is 1.
	takes as input a list of scores and returns a new list of scaled scores.
	'''
	# Is there really not a more elegant way to do this
	oldMin = min(scores)
	oldMax = max(scores)
	newMin = 1
	newMax = 10**selectivePressure
	newScores =  [(newMax - newMin)*(x - oldMin)/(oldMax - oldMin) + newMin for x in scores]
	return [s/sum(newScores) for s in newScores]

def optimizeMatrix_geneticAlg(scorematrix, mutationChance, mutationAmount,
	selectivePressure, N, totalStepsToStop, stepsWithNoImprovement,
	librarySize, gap_start, gap_extend, true_pos_align, true_neg_align):
	'''
	optimize an alignment score matrix using a genetic algorithm.
	stops if we hit totalItersToStop iterations, or if we don't see
	a new objective function value in the top (librarySize) in
	stepsWithNoImprovement steps.
	takes as input:
	scorematrix
	mutationChance		probability each entry in a scoring matrix will mutate
	mutationAmount		stdev of mutation gaussian
	selectivePressure	most fit individual is scaled to 10^selectivePressure * least fit
	n                   pop size
	totalStepsToStop	runtime step cutoff
	stepsWithNoImprovement
						step cutoff for seeing no new entries in our library of best matrices
	librarySize			number of best matrices to remember and re-seed population with
	gap_start			alignment params
	gap_extend
	true_pos_align		true positives
	true_neg_align		true negatives
	returns the final population of matrices, the scores for those matrices at the last step,
	the library of best matrices, and a list of the mean objective function value at each step.
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
