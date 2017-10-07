from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw,MCS
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Fingerprints import FingerprintMols
from sklearn.metrics.pairwise import pairwise_distances
import scipy
import scipy.cluster.hierarchy as sch  
from scipy.cluster import hierarchy
from numpy import arange
import numpy as np
import pandas as pd  
import matplotlib.pylab as plt
import os

# Note: MCSFP and MCFFP developed for the first time to calculating the Max Commen Substructures Fingerprints and Max Commen Frameworks Fingerprints
# Author information: shenwanxiang, any bugs are welcomed and users can send Emails to senwanxiang@163.com


# calculating Morgan FP
def CalMorganFp(radius,mols):
    # calculating the fingerprint
    MorganFP = []
    for mol in mols:
        fp = AllChem.GetMorganFingerprint(mol,radius)
        MorganFP.append(fp)
    return MorganFP

def CalTopologicalFp(mols):
    # calculating the fingerprint
    MorganFP = []
    for mol in mols:
        fp = FingerprintMols.FingerprintMol(mol)
        MorganFP.append(fp)
    return MorganFP

# calculating tanimoto similarity,shape:n*(n-1)/2
def TanimotoList(MorganFP):    
    Tanimoto_list =[]   
    for i in range(len(MorganFP)-1):
        for fp2 in MorganFP[i+1:]:
            s = DataStructs.TanimotoSimilarity(MorganFP[i],fp2) #Available similarity metrics include Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, and Tversky
            Tanimoto_list.append(s)
    return Tanimoto_list
 

# convert Z to newick formats
def GetNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = GetNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = GetNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
    return newick



def SmilesToClass(data):
    Smiles = list(data['smiles'])
    # smiles_list to mols
    mols = [Chem.MolFromSmiles(x) for x in Smiles]
    # calculating MorganFP and Tanimoto similarity
    MorganFP = CalMorganFp(4,mols)
    Tanimoto_list = TanimotoList(MorganFP)
    # convert to distance D
    D = 1-np.array(Tanimoto_list)
    # using the tanimoto distance to cluster
    Z= sch.linkage(D,method='complete')
    # convert dedrogram
    fig = plt.figure(figsize=(20,13))
    P =sch.dendrogram(Z,color_threshold=0.9,leaf_font_size=9,distance_sort = 'ascending') #labels=leaf_names, truncate_mode='lastp',p=12
    #plt.show()
    plt.savefig(str(0.9)+'.png',dpi = 300)
    plt.close()
    # Divided into diffrent classes, the threshold is set as Similarity_threshold
    return Z
    









if __name__ == '__main__':

    # set parameters
    start = 0.3
    stop = 0.95
    step=0.05
    data = pd.read_csv('data/scaffolds.txt')

    class_list = []
    for i in arange(start,stop,step):
        df_class = SmilesToClass(data,distance_threshold = i)
        print(i)
        class_list.append(df_class)
    #all_class = pd.concat(class_list.axis = 1)

'''
import pandas as pd
import SmilesToClass as STC
import numpy as np
from numpy import arange
import scipy.cluster.hierarchy as sch 

start = 0.3
stop = 0.95
step=0.05
data = pd.read_csv('data/scaffolds.txt')
Z = STC.SmilesToClass(data)
df = pd.DataFrame()
for i in arange(start,stop,step):
    T = sch.fcluster(Z, i, 'distance')
    # get the class merged into the dataframe
    df_class = pd.DataFrame({'group'+str(i):T})    
    print(i)
    df = pd.concat([df,df_class],axis=1)



    
    
lenght=[]
for i in df.columns:
    lengtht.append(len(df[i].unique()))
x = arange(start,stop,step)



all_class = pd.concat(class_list.axis = 1)

'''





