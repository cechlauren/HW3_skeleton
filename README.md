# Smith-Waterman

[![Build
Status](https://travis-ci.org/cechlauren/HW3_skeleton.svg?branch=master)](https://travis-ci.org/cechlauren/HW3_skeleton)

Smith-Waterman-derived alignment project with testing.

## usage

To use the package, first make a new conda environment and activate it

```
conda create -n exampleenv python=3
source activate exampleenv
```

then run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `HW3_skeleton/hw3align/__main__.py`) can be run as follows

```
python -m smith_waterman
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.

## Questions Part 1

### Question 1
Consider the false positive rate (proportion of negative pairs with scores that exceed
a score threshold) when the true positive rate (proportion of positive pairs with scores
above the threshold) is 0.7. 
-What's the best false positive rate that you can achieve with varying both gap opening (from 1 to 20) and extension penalties (from 1 to 5) with the BLOSUM50 matrix? 
-What is the best gap penalty combination?

The best (read lowest) false positive rate I could achieve was 20% when varying gap opening/extension given a true positive rate of 70%.  This is shown in [optimalgaps.py](https://github.com/cechlauren/HW3_skeleton/blob/master/hw3align/optimalgaps.py).
R was used to plot this data (a table of values fpr values for each combinatino of gap opening/extension), which can be found in [gapPenalties.txt](https://github.com/cechlauren/HW3_skeleton/blob/master/hw3align/gapPenalities.txt), to find the optimal combination of gap penalties. 
Here is the R script: [OptimizeGapsRscript.txt](https://github.com/cechlauren/HW3_skeleton/blob/master/OptimizeGapsRscript.txt).

We determine all true positive alignment scores to find the cutoff that sets true positive rate at 70%; they range from roughly 35-250.
The FPR will then be the count of neg. scores above the cutoff we designated, divided by the number of true negatives. This is all described in [optimalgaps.py](https://github.com/cechlauren/HW3_skeleton/blob/master/hw3align/optimalgaps.py).

So, using those cutoffs, we get the following false positive rate distribution:
<img src="falseposrate.png" /><br />
If you see nothing, also see: [optimizeGapPenalitiesPlot.pdf](https://github.com/cechlauren/HW3_skeleton/blob/master/optimizeGapPenalitiesPlot.pdf)

I found several combinations [lowestFPR.png](https://github.com/cechlauren/HW3_skeleton/blob/master/lowestFPR.png) that result in a 20% FPR, but decided to go with the less extreme case where

Gap opening: -8
Gap extension: -3


The following analyses will use the above affine penalties.

### Question 2
-Using the gap penalties you determined from question 1, which of the provided
scoring matrices performs the best, in terms of false positive rate (at a true positive rate
of 0.7)? 
-What are the performance rates of each of the matrices? 

To identify the best scoring matrix based on the best affine penalities I identified in question 1 I calculated and plotted the respective scores using [roc.py](https://github.com/cechlauren/HW3_skeleton/blob/master/hw3align/roc.py). 
All plots can be viewed in this directory: [ROCplots](https://github.com/cechlauren/HW3_skeleton/tree/master/ROCplots).
The following are the unnormalized ROC for each scoring matrix:
<img src="blosum50_8_3.png" /><br />

<img src="blosum62_8_3.png" /><br />

<img src="pam100_8_3.png" /><br />

<img src="pam50_8_3.png" /><br />

<img src="matio_8_3.png" /><br />






