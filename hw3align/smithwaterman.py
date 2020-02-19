import sys

import numpy as np




def sw(a, b, substitutionmatrix, startcost, extendcost):

	'''

	uses smith-waterman to find a maximal local alignment.
	we'll have a big table that has indices keeping track of the aligned 
	seq's beginning and ends. An ideal solution will only need to compute
	each sub-alignment once, while providing coverage of the whole
	search space.


	The input for this function is two strings (a,b) to align,  a

	substitution matrix, a gap start cost, and finally a gap extension cost.

	Most commonly in evolution, gaps will occur in bunches, so it would make 
	sense that once a break in alignment is made that it would be easier to 
	make multiple insertions or deletions. In contrast, to introduce the first 
	gap (the "start"/"open" gap) would be more costly since a break must occur
	in the DNA. 
	Therefore, a gap cost is: (d + (n-1)*e), where "n" is the length of

	the gap, "d" is the gap opening, and "e" is the gap extend.
	By this logic, the gap extension cost will only be only bepaid for 2+ 
	length gaps.



	This function will return:
	1) a tuple of the indices for beginning of local alignment in a and b 
	2) a tuple of the indices of the end of the alignment 
	3) the score of the aligning section 
	4) the aligned sections of a and b (gaps = "-").

	'''

	# initialize our matrices and trace matrices
	#Local alignment initializes with top left =0

	match = np.zeros(shape=(len(a) + 1, len(b) + 1))

	traceMatch = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)]

	Y = np.zeros(shape=(len(a) + 1, len(b) + 1))

	traceY = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)]

	X = np.zeros(shape=(len(a) + 1, len(b) + 1))

	traceX = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)]

	# we fill the first rows/columns of matrices to avoid alignments causing unusual edges

	for i in range(len(a) + 1):

		Y[i, 0] = float("-inf")

		if i>0: 

			X[i, 0] = startcost + extendcost*(i-1)

	for j in range(len(b) + 1):

		X[0, j] = float("-inf")

		if j>0:

			Y[0, j] = startcost + extendcost*(j-1)

	# fill in matrices

	maxSeen = (float("-inf"), (-1, -1)) # (score, (i,j))

	for i in range(1, len(a) + 1):

		for j in range(1, len(b) + 1): 

			extendX = X[i-1, j] + extendcost  # extend gap in b

			openX = match[i-1, j] + startcost # open gap in b

			if openX > extendX:

				X[i,j] = openX

				traceX[i][j] = "openX"

			else:

				X[i,j] = extendX

				traceX[i][j] = "extendX"

			extendY = Y[i, j-1] + extendcost  # extend gap in a

			openY = match[i, j-1] + startcost # open gap in a

			if openY > extendY: 

				Y[i,j] = openY

				traceY[i][j] = "openY"

			else:

				Y[i,j] = extendY

				traceY[i][j] = "extendY"

			match[i,j] = match[i-1, j-1] + substitutionmatrix[frozenset((a[i-1], b[j-1]))]

			traceMatch[i][j] = "match"

			if Y[i,j] > match[i,j]:

				match[i,j] = Y[i,j]

				traceMatch[i][j] = "closeY"

			if X[i,j] > match[i,j]:

				match[i,j] = X[i,j]

				traceMatch[i][j] = "closeX"

			if match[i,j] > maxSeen[0]:

				maxSeen = (match[i,j], (i, j))

	# traceback

	i,j = maxSeen[1]

	end = (i,j)

	matrix, traceMatrix = match, traceMatch

	score = matrix[i,j]

	bestScore = score

	matchedA = ""

	matchedB = ""

	while score > 0:

		if traceMatrix[i][j] == "match":

			matchedA = a[i-1] + matchedA

			matchedB = b[j-1] + matchedB

			i,j = i-1, j-1

		elif traceMatrix[i][j] == "closeX":

			matrix, traceMatrix = X, traceX

		elif traceMatrix[i][j] == "closeY":

			matrix, traceMatrix = Y, traceY

		elif traceMatrix[i][j] == "openY":

			matchedA = "-" + matchedA

			matchedB = b[j-1] + matchedB

			i,j = i, j-1

			matrix, traceMatrix = match, traceMatch

		elif traceMatrix[i][j] == "extendY":

			matchedA = "-" + matchedA

			matchedB = b[j-1] + matchedB

			i,j = i, j-1

		elif traceMatrix[i][j] == "openX":

			matchedA = a[i-1] + matchedA

			matchedB = "-" + matchedB

			i,j = i-1, j

			matrix, traceMatrix = match, traceMatch

		elif traceMatrix[i][j] == "extendX":

			matchedA = a[i-1] + matchedA

			matchedB = "-" + matchedB

			i,j = i-1, j

		else: print("we should never get here.")

		score = matrix[i,j]

	start = (i, j)

	#print(start, end, bestScore)

	#print(matchedA)

	#print(matchedB)

	return start, end, bestScore, matchedA, matchedB




