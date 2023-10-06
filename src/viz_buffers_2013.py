'''
viz_buffers_2013.py

Trying to match HOOH production data from Morris et al 2011 fig 1 using an analytical solution, euler's aproximation, and ODEint 

UH18301 Pro in 3.75 mM HEPES or Taps buffer. High light (24uC in a Sunbox -  noon maximum of about 250 quanta m) pro99 media 

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




###############################

#  Importing the data frame

###############################

 #UH18301 Pro in 3.75 mM HEPES or Taps buffer. High light (24uC in a Sunbox -  noon maximum of about 250 quanta m) 


df_all = pd.read_csv('../data/all_buffers_2013Morrisf1c.csv')

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#creating log and stats for data 
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])
df_all['log4'] = np.log(df_all['rep4'])
df_all['log5'] = np.log(df_all['rep5'])
df_all['log6'] = np.log(df_all['rep6'])


df_all['avg1'] = np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['avg2'] = np.nanmean(np.r_[[df_all[i] for i in ['rep4', 'rep5', 'rep6']]],axis=0)
df_all['std1'] = np.std(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['std2'] = np.std(np.r_[[df_all[i] for i in ['rep4', 'rep5', 'rep6']]],axis=0)

#bio rep logged avgs and stdvs 
df_all['lavg1'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['lavg2'] = np.nanmean(np.r_[[df_all[i] for i in ['log4','log5','log6']]],axis=0)
df_all['stdlog1'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['stdlog2'] = np.std(np.r_[[df_all[i] for i in ['log4','log5','log6']]],axis=0)



df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4', 'rep5', 'rep6']]],axis=0)
df_all['sigma'] =  np.std(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4', 'rep5', 'rep6']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3','log4','log5','log6']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3','log4','log5','log6']]],axis=0)


##############################

#   data arrays

###############################


orgs = df_all['organism'].unique()

assays = df_all['ID'].unique()
nassays = assays.shape[0]
treats  = df_all['Seawater'].unique()
ntreats = treats.shape[0]

#HOOHs_nm = 1000**np.array(HOOH_data)   #TRIED TO GET IN nM to match better with other things 
#print(HOOH_data)
colors = ('r','g')
markers = ('o','*')

##############################

#    graphing the data 

##############################


for a,n in zip(assays,range(nassays)):
    fig1,(ax1)= plt.subplots(2,1, figsize = (12,8))
    fig2,(ax2) = plt.subplots(1,2,figsize = (12,8))
    for t,nt in zip(treats,range(ntreats)): 
        df = df_all[(df_all['ID'] == a)]
        count = nt
        df = df[(df['ID']==a) & (df['Seawater']==t)]
        pdf = df[df['organism'] == 'P']
        hdf = df[df['organism'] == 'H']
        ax1[0].plot(pdf['time'], pdf['abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        ax1[1].plot(hdf['time'], hdf['abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        fig1.suptitle('Abiotic Dynamics '+ str(a))
        ax1[1].set_ylabel('HOOH concentration (\u03BCM)')
        ax1[0].set_ylabel('Pro abundance (cell/ml)')
        fig1.supxlabel('Time (days)')
        l1 = ax1[0].legend(loc = 'lower right', prop={"size":14}) 
        l1.draw_frame(False)#print(df)
        ax1[0].semilogy()
        ax1[1].semilogy()
        #viz uncertainty in data
        ax2[0].errorbar(x = pdf['time'], y =pdf['log_abundance'],yerr = pdf['log_sigma'], color = colors[count], label = ('P in '+ str(t)))
        ax2[0].errorbar(x = hdf['time'], y =hdf['log_abundance'],yerr = hdf['log_sigma'], color = colors[count], label = ('H of '+ str(t)))
        ax2[1].plot(pdf['log_abundance'],pdf['log_sigma'], label = 'P')
        ax2[1].plot(hdf['log_abundance'],hdf['log_sigma'], label = 'H')
        fig2.suptitle('Mean and std of dynamics of '+ str(a))
        #suptitle('mean vs std')
        l2 = ax2[0].legend(loc = 'lower right', prop={"size":14}) 
        l2.draw_frame(False)#print(df)
        ax2[0].semilogy()
        ax2[0].set_ylabel('Concentration (ml-1)'+ str(a) +' in '+ str(t))
        ax2[0].set_xlabel('Time (Days)')
        ax2[1].set_xlabel('mean')
        ax2[1].set_ylabel('std')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show()
    fig1.savefig('../figures/Hproduction_'+str(a)+'_data.png')
    fig2.savefig('../figures/dynamics_'+str(a)+'_data.png')










    
'''


'''
print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
