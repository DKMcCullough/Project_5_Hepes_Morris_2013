'''
viz_seawater_hepes.py

Hepes 0 or 10 in different concentrations of seawater to miliQ water. 3 tech reps. 65 light setting

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

#H data only 

df_all = pd.read_csv('../data/seawater_assay_2013.csv')

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_all.fillna(0)
#creating log and stats for data 
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

#making avg and stdv of 3 reps 

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['sigma'] =  np.nanstd(np.r_[[df_all[i] for i in ['rep1','rep2','rep3',]]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.nanstd(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)


##############################

#   data arrays

###############################


assays = df_all['Buffer'].unique()
nassays = assays.shape[0]

treats  = df_all['Seawater.Percent'].unique()
ntreats = treats.shape[0]


##############################

#    graphing the data 

##############################


for a,n in zip(assays,range(nassays)):
    dfw = df_all[(df_all['Buffer'] == a)]
    #print(dfw)
    fig1,(ax1)= plt.subplots(ntreats,2, figsize = (10,8))
    fig1.suptitle('H raw Dynamics from'+ str(a) +'% Seawater', size = 22)
    fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
    ax1[2,0].set_ylabel('HOOH concentration (\u03BCM)')
    ax1[4,0].set_xlabel('Time' )
    ax1[0,0].semilogy()
    ax1[2,1].set_ylabel('STDV')
    ax1[4,1].set_xlabel('Mean' )
    fig2,(ax2) = plt.subplots(ntreats,2,figsize = (10,8))
    fig2.suptitle(' Log  dynamics of '+ str(a)+' % Seawater', size = 22)
    fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
    ax2[2,0].set_ylabel('HOOH concentration (\u03BCM)')
    ax2[4,0].set_xlabel('Time' )
    ax2[0,0].semilogy()
    ax2[2,1].set_ylabel('Log STDV')
    ax2[4,1].set_xlabel('Log Mean' )
    for t,nt in zip(treats,range(ntreats)):
        df = dfw[(dfw['Seawater.Percent']==t) & (dfw['Buffer']==a)]
        #print off reps and then raw mean and stdv
        ax1[nt,0].plot(df['time'], df['rep1'], marker= '*', markersize= 10, label =('rep1'), color = 'pink') 
        ax1[nt,0].plot(df['time'], df['rep2'], marker= 'o', markersize= 10, label =('rep2'), color = 'purple' ) 
        ax1[nt,0].plot(df['time'], df['rep3'], marker= 'd', markersize= 10, label =('rep3'), color = 'c' ) 
        ax1[nt,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker= '.',markersize= 10, label ='Raw Mean', color = 'red' )
        ax1[nt,0].text(1.2,0.5,'Seawater % : '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
        #give dynamic graph rep legend and log it
        ax1[nt,0].semilogy()
        l1 = ax1[nt,0].legend(loc = 'lower center')
        l1.draw_frame(False)
        #stdv vs mean of raw data 
        ax1[nt,1].scatter(df.abundance,df.sigma, c='purple')
        #print out logged data and logged mean
        ax2[nt,0].plot(df['time'], df['log1'], marker= '*', markersize= 10, label =('log1'), color = 'pink') 
        ax2[nt,0].plot(df['time'], df['log2'], marker= 'o', markersize= 10, label =('log2'), color = 'purple' ) 
        ax1[nt,0].plot(df['time'], df['log3'], marker= 'd', markersize= 10, label =('log3'), color = 'c' ) 
        ax2[nt,0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma, marker = '*', c='red',label =  'Log Mean')
        #log y axis and give rep legend         
        l2 = ax2[nt,0].legend(loc = 'lower right')
        l2.draw_frame(False)
        ax2[nt,0].semilogy()
        
        #annotate side of graph with treatment 
        ax2[nt,0].text(1.2,0.5,'Seawater % : '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
#give bigger graph ticks and print out graphs 
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show()
    #save figures 
    fig1.savefig('../figures/Hproduction_'+str(a)+'seawater_data.png')
    fig2.savefig('../figures/dynamics_'+str(a)+'seawater_data.png')










    
'''


'''
print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
