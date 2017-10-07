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


def GetMCSSmiles(mol,labelledMol,scs):
    # input mol, copyed mol, smart string
    # output smiles of the Maximum Common Substructure
    mcsp = Chem.MolFromSmarts(scs)
    match = labelledMol.GetSubstructMatch(mcsp)
    return Chem.MolFragmentToSmiles(mol,atomsToUse=match,
                                    isomericSmiles=False,
                                    canonical=False)

def GetMCSNum(df, class_group = 'class', smiles_type = 'smiles'):
    # input: Dataframe of T, smiles, scaffolds and frameworks
    # output: Dataframe of MCS and number of MCS              
    commen_substrc = []
    for i in set(df[class_group]):
        #print(i)
        smiles = df[df[class_group] == i]
        ms = [Chem.MolFromSmiles(x) for x in smiles[smiles_type]]
        num = len(ms)
        if num == 1:
            #Draw.MolToFile(ms[0],'result/class_'+str(i)+'_num_'+str(len(ms))+'.png')
            mcs = list(smiles[smiles_type])[0]
            scs = mcs
        else:
            #img=Draw.MolsToGridImage(ms,molsPerRow=4,subImgSize=(500,500))
            #img.save('result/class_'+str(i)+'_num_'+str(len(ms))+'.png')
            nms = [Chem.Mol(x) for x in ms] #copy ms since we will change them
            res = MCS.FindMCS(ms,timeout=10) #return the best match found in this time
            scs = res.smarts # smartstring 
            if scs == None:
                mcs = 'A'
            else: mcs = GetMCSSmiles(ms[1],nms[1],scs) #smiles
        cs_num=[mcs,scs,num,i]
        commen_substrc.append(cs_num)
    c_a = np.array(commen_substrc)
    df_result = pd.DataFrame({'MCS':c_a[:,0],'smarts':c_a[:,1],'num':c_a[:,2].astype('int64'),'class':c_a[:,3].astype('int64')})
    return df_result


def SmilesToClass(data,distance_threshold = 0.8,mcs = True):
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
    P =sch.dendrogram(Z,color_threshold=distance_threshold,leaf_font_size=9,distance_sort = 'ascending') #labels=leaf_names, truncate_mode='lastp',p=12
    #plt.show()
    try: os.mkdir('result')
    except OSError: pass
    plt.savefig('result/'+str(distance_threshold)+'.png',dpi = 300)
    plt.close()
    # Divided into diffrent classes, the threshold is set as Similarity_threshold
    T = sch.fcluster(Z, distance_threshold, 'distance')
    # get the class merged into the dataframe
    my_class = pd.DataFrame({'class':T})
    df_class = pd.concat([data,my_class],axis = 1) # concat good smiles and the new class
    #get the MCS of each class and draw the compunds
    if mcs == True:
        df_mcs = GetMCSNum(df_class, class_group = 'class', smiles_type = 'smiles')
        df_mcs.sort('num',ascending = False).to_csv('result/'+str(distance_threshold)+'.csv',index=None)
    else: df_mcs = ''
    return df_class,df_mcs



def McsToFp(smiles_list,mcs_list,fp_name_list):
    mols = [Chem.MolFromSmiles(x) for x in smiles_list]
    fp_name = np.array(fp_name_list).reshape(1,len(fp_name_list))
    m = len(mols)
    n = len(mcs_list)
    MY_FP = np.zeros([m,n])
    for i in range(len(mols)):
        #Draw.MolToFile(mols[i],'result/orignal/'+str(i+1)+'.png')
        for j in range(len(mcs_list)):
            print('this is mol %d and fp %d' % (i,j))
            try:
                if mols[i].HasSubstructMatch(Chem.MolFromSmiles(mcs_list[j],sanitize=False)):
                    MY_FP[i,j] = 1
            except: pass
    FP = np.concatenate((fp_name,MY_FP),axis = 0)
    return FP


def CleanFp(df_fp, df_mcs_id, distance_threshold = 0.1):
    # this func aims to cluster the MCS or MCF fingerprints by many moleculars
    # imput: DataFrame of fingerprints for many mols :df_fp, df for mcs and ids with no head: df_mcs_id
    # out_put: clearned fingerprints, the duplicated fps are removed and smiliar fps are convert to the minist ones
    fp_names = list(df_fp.columns)
    fp_values = df_fp.values.astype('int64').T    
    # calculating tanimoto distance list 
    D_list = scipy.spatial.distance.pdist(fp_values,'rogerstanimoto') # list with length n*(n-1)/2, VIP * X = squareform(v) list convert to matrix or reverse* 'correlation'
    # using the tanimoto distance to cluster
    Z= sch.linkage(D_list,method='complete')
    # convert to newick format tree
    tree = hierarchy.to_tree(Z,False)
    tree.newick = GetNewick(tree, "", tree.dist, fp_names)
    with open('tree.newick','w') as f:
        f.write(tree.newick)
    # convert dedrogram
    fig = plt.figure(figsize=(15,10))
    P =sch.dendrogram(Z,color_threshold=distance_threshold,leaf_font_size=6,labels =fp_names,distance_sort = 'ascending') #labels=leaf_names, truncate_mode='lastp',p=12
    plt.show()
    #plt.savefig(str(distance_threshold)+'_fp.png',dpi = 100)
    plt.close()
    # Divided into diffrent classes, the threshold is set as Similarity_threshold
    T = sch.fcluster(Z, distance_threshold, 'distance')
    # get the class merged into the dataframe
    df_class = pd.DataFrame({'class':T,'id':fp_names})
    gg = pd.concat([df_mcs_id,df_class],axis = 1)
    id_list = []
    for class_ in set(gg[gg.columns[2]]):
        oneclass = gg[gg['class'] == class_]
        if len(oneclass) == 1:
            oneid = list(oneclass['id'])[0]
        else:
            t = oneclass[oneclass.columns[0]]
            kk = [len(x) for x in t].index(min([len(x) for x in t]))
            ss = list(oneclass['id'])
            oneid = ss[kk]
        id_list.append(oneid)
    df_select = df_mcs_id[df_mcs_id[df_mcs_id.columns[1]].isin(id_list)]
    return df_class,df_select





if __name__ == '__main__':
    # read the mother template
    smiles_list = list(pd.read_csv('template.txt',sep='\t',header =None)[0])

    # read the self-defined fps
    fp = pd.read_csv('mcs.smiles',sep='\t',header =None)
    mcs_list = list(fp[0])
    fp_name_list = list(fp[1])

    #calculation the fingerprints
    FFP = McsToFp(smiles_list,mcs_list,fp_name_list)
    np.savetxt('calculated_fp.csv',FFP,fmt='%s',delimiter=',')










