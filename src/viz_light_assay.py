'''
viz_light_assay.py


Visualizing the light assay from Morris 2013 (fig 1c) . 0 or 10M hepes in differing amouuns of light. !00% seawater 

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


df_all = pd.read_csv('../data/light_assay_2013.csv')

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_all.fillna(0)
#creating log and stats for data 
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

#making avg and stdv of 3 reps 

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['sigma'] =  np.nanstd(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.nanstd(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)

##############################

#   data arrays

###############################



assays = df_all['Buffer'].unique()
nassays = assays.shape[0]

treats = df_all['Light'].unique()
ntreats = treats.shape[0]



##############################

#    graphing the data 

##############################

for a,n in zip(assays,range(nassays)):
    dfw = df_all[(df_all['Buffer'] == a)]
    fig1,(ax1)= plt.subplots(ntreats,2, figsize = (14,11))
    fig1.suptitle('H raw Dynamics from'+ str(a) +' light', size = 22)
    fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
    ax1[3,0].set_ylabel('HOOH concentration (\u03BCM)', size = 16)
    ax1[5,0].set_xlabel('Time' , size = 16)
    ax1[0,0].semilogy()
    ax1[3,1].set_ylabel('STDV', size = 16)
    ax1[5,1].set_xlabel('Mean' , size = 16)
    fig2,(ax2) = plt.subplots(ntreats,2,figsize = (14,11))
    fig2.suptitle(' Log  dynamics of '+ str(a), size = 22)
    fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
    ax2[3,0].set_ylabel('HOOH concentration (\u03BCM)', size = 16)
    ax2[5,0].set_xlabel('Time' , size = 16)
    ax2[0,0].semilogy()
    ax2[3,1].set_ylabel('Log STDV', size = 16)
    ax2[5,1].set_xlabel('Log Mean', size = 16 )
    for t,nt in zip(treats,range(ntreats)):
        df = dfw[(dfw['Light']==t) & (dfw['Buffer']==a)]
        ax1[nt,0].plot(df['time'], df['rep1'], marker= '*', markersize= 10, label =('rep1'), color = 'pink') 
        ax1[nt,0].plot(df['time'], df['rep2'], marker= 'o', markersize= 10, label =('rep2'), color = 'purple' ) 
        ax1[nt,0].plot(df['time'], df['rep3'], marker= 'd', markersize= 10, label =('rep3'), color = 'c' ) 
        ax1[nt,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker= '.',markersize= 10, label =(str(t)+' Mean produced HOOH'), color = 'brown' )
        ax1[nt,0].text(1.2,0.5,'Light: '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
        ax1[nt,0].semilogy()
        l1 = ax1[nt,0].legend(loc = 'upper left')
        l1.draw_frame(False)
        ax1[nt,1].scatter(df.abundance,df.sigma, c='purple')
        ax2[nt,0].semilogy()
        ax2[nt,0].plot(df['time'], df['log1'], marker= '*', markersize= 10, label =('log1'), color = 'pink') 
        ax2[nt,0].plot(df['time'], df['log2'], marker= 'o', markersize= 10, label =('log2'), color = 'purple' ) 
        ax1[nt,0].plot(df['time'], df['log3'], marker= 'd', markersize= 10, label =('log3'), color = 'c' ) 
        ax2[nt,0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma, marker = '*', c='brown',label =  'Log Mean produced HOOH')
        #annotate
        ax2[nt,0].text(1.2,0.5,'Light: '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
        l2 = ax2[nt,0].legend(loc = 'upper left')
        l2.draw_frame(False)
        #stdv vs mean


plt.show()
fig1.savefig('../figures/Hproduction_'+str(a)+'_light.png')
fig2.savefig('../figures/dynamics_'+str(a)+'_data.png')






    
'''


'''
print('\n ~~~****~~~****~~~ \n')
print('done with light assay ')
print('\n ~~~****~~~****~~~ \n')
