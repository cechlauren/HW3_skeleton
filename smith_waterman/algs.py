import numpy as np
import pandas as pd
#from .utils import Atom, Residue, ActiveSite
#from .io import read_active_sites, read_active_site, write_clustering, write_mult_clusterings

"""
###########################################################################################################
#PROVIDE CONTEXT TO MY METHODS:
# I consider each residue in my active sites to embody 4 different physiochemical properties--
#this means each of the 20 amino acids is now condensed to an alphabet of just 6 letters.
#these properties are: L = large, A= acidic, B= basic, N= neutral, F= aliphatic, P= polar
#By reducing the features of my active sites from 20 to 6, I hope to find trends in composition across active sites.
#Similar composition and similar size would suggest similar ligands (the things that bind to active sites), and thus 
#similar biological function.
###########################################################################################################
#MAKE A DICTIONARY:
#Here's a dictionary of those residues with their corresponding physiochemical definitions:
convert = {'ARG': 'B', 'LYS': 'B', 'ASP': 'A', 'GLN': 'P', 'ASN': 'P', 'GLU': 'A',
'HIS': 'B', 'SER': 'P', 'THR': 'P', 'PRO': 'N', 'TYR': 'P', 'CYS': 'P', 'GLY': 'N', 'ALA': 'F',
'MET': 'N', 'TRP': 'L', 'LEU': 'F', 'VAL': 'A', 'PHE': 'L', 'ILE': 'F'}
###########################################################################################################
#PREP MY ACTIVE SITE DATA:
#Now I need to extract all of those active sites from the data folder, using those utility functions.
#First call the function that reads all active sites in a folder.
#I extract the residues from each active site, make them a list, then add them to some temporary list of all active sites.
#Next I need to change each of those residue names into my own kind.
#The residues are still separated, and I need them in one string, so join them together.
#Finally, I get them into the format I used to start making my similarity matrix. 
My_active_Site = read_active_sites("../data/")
temp = []
for site in My_active_Site:
    temp_physio = []
    amino_acid_list = site.residues
    for amino_acid in amino_acid_list:
        temp_physio.append(amino_acid.type)
    temp.append(temp_physio)
new_temp = []
for item in temp:
    new_temp_physio = []
    for aa in item:
        new_temp_physio.append(convert[aa])
    new_temp.append(new_temp_physio)
new_super_temp = []
for almost in new_temp:
    new_super_temp.append(''.join(almost))
my_res_list = [ [x] for x in new_super_temp]
# print(my_lame_list) #if you want
#now I should have a list of all my active sites' residues, composed of only those 6 specified features.
###########################################################################################################
#MAKE A DATAFRAME OF MY ACTIVESITES X FEATURES:
#now I need to make my dataframe with all my final values of the counts of those residues in my active sites.
#(Miriam helped with this):
#Using map allows us to apply a function to directly to a list.
counter_list=list(map(lambda x: list(map(lambda y: Counter(y),x)),my_res_list)) #This is a counter for each of the features
list_of_dfs=list(map(lambda z: pd.DataFrame.from_dict(z),counter_list)) #I think this makes a data frame for each of the counts
dfObj = pd.DataFrame() #make an empty data frame to add our count data frame to
a_df = dfObj.append(list_of_dfs, sort = False).fillna(0) #join those two data frames together and fill in all NAs with zero
#This outputs a 136xfeature dataframe with each feature having some number of counts based on the activesite. 
#My activesites are not labeled, so need to indexed them
#a_df.index=range(len(my_res_list)) #index them based on length of the list of active sites I put in
active_site_index = [46495,23812,41729,91911,82212,15813,85232,20856,3458,13052,82993,32088,64392,29047,42633,53272,97218,69893,96099,82238,50018,68578,55996,93456,33838,71389,43878,57602,58445,39939,81697,17622,81859,6040,57370,22711,52235,46042,34563,78796,9776,83741,7674,91796,63064,24307,63703,8208,62186,28919,56394,81816,34088,37224,40084,82886,23319,42296,93168,42269,17526,38472,83227,29209,83394,10701,94372,93192,98170,85492,73462,38846,91426,27031,28672,46975,25551,91194,18773,37237,25196,61242,56029,42074,70005,35014,19267,1806,57644,26095,64258,84035,45127,34047,14181,57481,29773,54203,36257,4629,26246,37438,98797,81563,65815,38181,63634,23760,49624,39117,24634,94652,7780,73624,3733,73183,42202,32054,50362,276,70919,94719,10814,25878,39299,27312,88042,38031,52954,20326,8304,72058,34958,41719,97612,47023]
a_df.index= active_site_index #I now have indexed by df using the correct active site names
#I also included a total count list that can inform me about the size of the active site itself
a_df["total"]=a_df.sum(axis=1)
#Determine the similarity between all given ActiveSite instances, where the input is all active site's features
similarity_df = a_df.T.corr()
#(residues and total) counts, and the output will be the "distance" or dissimilarity between them (since I do 1-correlation).
distance_matrix = 1-a_df.T.corr()
###########################################################################################################
"""


#I'll consider my similarity metric to be the "correlation" between active site instances that I produced above.
#I'll copy and paste all of those matrix construction machinations here for the similarity function.
"""
def compute_similarity(site_a, site_b):
  
    
    convert = {'ARG': 'B', 'LYS': 'B', 'ASP': 'A', 'GLN': 'P', 'ASN': 'P', 'GLU': 'A', 'HIS': 'B', 'SER': 'P', 'THR': 'P', 'PRO': 'N', 'TYR': 'P', 'CYS': 'P', 'GLY': 'N', 'ALA': 'F', 'MET': 'N', 'TRP': 'L', 'LEU': 'F', 'VAL': 'A', 'PHE': 'L', 'ILE': 'F'}
    
    My_active_Site = read_active_sites("../data/")
    temp = []
    for site in My_active_Site:
        temp_physio = []
        amino_acid_list = site.residues
        for amino_acid in amino_acid_list:
            temp_physio.append(amino_acid.type)
        temp.append(temp_physio)
    new_temp = []
    for item in temp:
        new_temp_physio = []
        for aa in item:
            new_temp_physio.append(convert[aa])
        new_temp.append(new_temp_physio)
    new_super_temp = []
    for almost in new_temp:
        new_super_temp.append(''.join(almost))
    my_res_list = [ [x] for x in new_super_temp]
    
    counter_list=list(map(lambda x: list(map(lambda y: Counter(y),x)),my_res_list))
    list_of_dfs=list(map(lambda z: pd.DataFrame.from_dict(z),counter_list))
    dfObj = pd.DataFrame()
    a_df = dfObj.append(list_of_dfs).fillna(0)
    
    active_site_index = [46495,23812,41729,91911,82212,15813,85232,20856,3458,13052,82993,32088,64392,29047,42633,53272,97218,69893,96099,82238,50018,68578,55996,93456,33838,71389,43878,57602,58445,39939,81697,17622,81859,6040,57370,22711,52235,46042,34563,78796,9776,83741,7674,91796,63064,24307,63703,8208,62186,28919,56394,81816,34088,37224,40084,82886,23319,42296,93168,42269,17526,38472,83227,29209,83394,10701,94372,93192,98170,85492,73462,38846,91426,27031,28672,46975,25551,91194,18773,37237,25196,61242,56029,42074,70005,35014,19267,1806,57644,26095,64258,84035,45127,34047,14181,57481,29773,54203,36257,4629,26246,37438,98797,81563,65815,38181,63634,23760,49624,39117,24634,94652,7780,73624,3733,73183,42202,32054,50362,276,70919,94719,10814,25878,39299,27312,88042,38031,52954,20326,8304,72058,34958,41719,97612,47023]
    a_df['']= active_site_index
    b_df = a_df.set_index('')
    
    
    b_df["total"]=b_df.sum(axis=1) #cant use a_df sum bc will include index values...
    
    similarity_df = b_df.T.corr()
    
    
    #now I have a matrix containing all of the correlations between activesites. 
    #The higher the score (0:1), the higher the similarity between them. 
    #To get the similarity, find the instance in the similarity matrix when the two activesites intersect.
    similarity = similarity_df[site_a][site_b]
    return similarity





"""

def readsequence():
    return None

def sw():

    return None


def score():
    return None

def roc():
    return None
