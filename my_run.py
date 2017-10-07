import pandas as pd
import numpy as np
import MCSF_FP as MCSFP
from numpy import arange



############################# Hunt the MCS or MCF ###################################

# set parameters
start = 0.3
stop = 0.95
step=0.05
min_num = 2 #each class to find MCSFP must have at least 2 drugs
data = pd.read_csv('data/scaffolds_test.txt',sep='\t')

mcs_dict = {}
for i in arange(start,stop,step):
    if i == start:
        df_class,df_mcs = MCSFP.SmilesToClass(data,distance_threshold = i)
        mcs_dict[i] = df_mcs
        print(i)     
    else:
        retrive_class = list(df_mcs[df_mcs['num'] < min_num]['class'])
        df_retrive = df_class[df_class['class'].isin(retrive_class)]
        data_new = df_retrive.rename(columns = {'class':'old_class'+str(i-step)}).reset_index(drop=True)
        df_class,df_mcs = MCSFP.SmilesToClass(data_new,distance_threshold = i)
        mcs_dict[i] = df_mcs
        print(i)

Dis = arange(start,stop,step)
mcs_fp = pd.concat([mcs_dict[x] for x in Dis],axis = 0)
mcs_fp[mcs_fp['num'] > min_num].to_csv('result/mcs_gt'+str(min_num)+'.csv',index=None)




############################# deal with MCS_MCF FP ###################################

# read the mother template
smiles_list = list(pd.read_csv('data/frameworks.txt',sep='\t',header =None)[0])

# read the self-defined fps
fp = pd.read_csv('data/MCF2FP.smiles',sep='\t')
mcs_list = list(fp[fp.columns[0]])
fp_name_list = list(fp[fp.columns[1]])

# calculation the fingerprints
FFP = MCSFP.McsToFp(smiles_list,mcs_list,fp_name_list)
np.savetxt('MCF2FP_calculated_fp.csv',FFP,fmt='%s',delimiter=',')


# clean the fingerprints
df = pd.read_csv('MCF2FP_calculated_fp.csv')
ss = df.sum()
ss.sort()
fp_num = ss.to_frame(name='times')

re_list = list(fp_num[fp_num['times']==0].index)

df_class, df_select = MCSFP.CleanFp(df, fp, distance_threshold = 0.002)
df_class.to_csv('MCF2FP_class.csv',index=None)
df_select.to_csv('MCF2FP_select.csv',index=None)
############################# deal with MCS_MCF FP ###################################

print('good boy!')




