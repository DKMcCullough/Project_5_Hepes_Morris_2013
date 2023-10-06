'''
viz_miliQ_buffers.py

visualise different buffer proiductions of H from Morris 2013 (miliQ_bufferH_2013MorrisSIfig1b) 
  HOOH production in 10 mM Buffers plus seawater media incubated  in lightexposed milli-Q water (i.e., without seawater solutes). rep1-3 are bio1 rep4-6 are bio2

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


df_all = pd.read_csv('../data/miliQ_bufferH_2013MorrisSIfig1b.csv')


df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'

df = df_all 

#logging data for latter graphing 
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])
df['log4'] = np.log(df['rep4'])
df['log5'] = np.log(df['rep5'])
df['log6'] = np.log(df['rep6'])
#####

#raw avgs reps

df['avg1'] =  np.nanmean(np.r_[[df[i] for i in ['rep1','rep2','rep3']]],axis=0)
df['avg2'] =  np.nanmean(np.r_[[df[i] for i in ['rep4','rep5', 'rep6']]],axis=0)
df['stdv1'] = np.nanstd(np.r_[[df[i] for i in ['rep1','rep2','rep3']]],axis=0)
df['stdv2'] = np.nanstd(np.r_[[df[i] for i in ['rep4','rep5', 'rep6']]],axis=0)
#total avgs and stdvs

df['abundance'] =  np.nanmean(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4','rep5', 'rep6']]],axis=0)
df['sigma'] = np.nanstd(np.r_[[df[i] for i in ['rep1','rep2','rep3','rep4','rep5', 'rep6']]],axis=0)

#log avgs of reps
df['lavg1'] = np.nanmean(np.r_[[df[i] for i in ['log1','log2','log3']]],axis=0)
df['lavg2'] = np.nanmean(np.r_[[df[i] for i in ['log4','log5', 'log6']]],axis=0)
df['stdlog1'] =  np.nanstd(np.r_[[df[i] for i in ['log1','log2','log3']]],axis=0)
df['stdlog2'] =  np.nanstd(np.r_[[df[i] for i in ['log4','log5', 'log6']]],axis=0)

#total avgs and stdvs

df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['log1','log2','log3','log4','log5', 'log6']]],axis=0)
df['log_sigma'] =  np.nanstd(np.r_[[df[i] for i in ['log1','log2','log3','log4','log5', 'log6']]],axis=0)


##############################

#   data arrays

###############################

treats = df_all['Buffer'].unique()
ntreats = treats.shape[0]


##############################

#    graphing the data 

##############################


fig1,(ax1)= plt.subplots(ntreats,2, figsize = (12,10))
fig1.suptitle('Raw H dynamics of Buffers in miliQwater')
fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax1[1,0].set_ylabel('HOOH concentration (\u03BCM)')
ax1[3,0].set_xlabel('Time' )
ax1[0,0].semilogy()
ax1[1,1].set_ylabel('STDV')
ax1[3,1].set_xlabel('Mean' )

fig2,(ax2) = plt.subplots(ntreats,2,figsize = (12,10))
fig2.suptitle(' Log  H dynamics of Buffers in miliQwater')
fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax2[1,0].set_ylabel('HOOH concentration (\u03BCM)')
ax2[2,0].set_xlabel('Time' )
ax2[0,0].semilogy()
ax2[1,1].set_ylabel('Log STDV')
ax2[3,1].set_xlabel('Log Mean' )

for t,nt in zip(treats,range(ntreats)):
    df = df_all[(df_all['Buffer'] == t)]
    ax1[nt,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker= '.',markersize= 10, label =('Mean'), color = 'c' )
    ax1[nt,0].errorbar(df.time,df.avg1, yerr=df.stdv1, marker= '.',markersize= 10, label =('avg1'), color = 'pink' )
    ax1[nt,0].errorbar(df.time,df.avg2, yerr=df.stdv2, marker= '.',markersize= 10, label =('avg2'), color = 'purple' )
    #annotate side of large fig
    ax1[nt,0].text(1.2,0.5,'Buffer: '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
    #log y axis and add legend to dynamics graph 
    ax1[nt,0].semilogy()
    l1 = ax1[nt,0].legend(loc = 'lower center')
    l1.draw_frame(False)
    #graph stats
    ax1[nt,1].scatter(df.abundance,df.sigma, c='c')
    ax1[nt,1].scatter(df.avg1,df.stdv1, c='pink')
    ax1[nt,1].scatter(df.avg2,df.stdv2, c='purple')
    #graph logged dynamics 
    ax2[nt,0].errorbar(df.time,df.log_abundance, yerr=df.sigma, marker= '.',markersize= 10, label =('Log Mean'), color = 'c' )
    ax2[nt,0].errorbar(df.time,df.lavg1, yerr=df.stdlog1, marker= '.',markersize= 10, label =('Log avg1'), color = 'pink' )
    ax2[nt,0].errorbar(df.time,df.lavg2, yerr=df.stdlog2, marker= '.',markersize= 10, label =('Log avg2'), color = 'purple' )
    #log y axis and add legend to dynamics graph 
    ax2[nt,0].semilogy()
    l2 = ax2[nt,0].legend(loc = 'upper right')
    l2.draw_frame(False)
    #graph loggedstats
    ax2[nt,1].scatter(df.log_abundance,df.log_sigma, c='c')
    ax2[nt,1].scatter(df.lavg1,df.stdlog1, c='pink')
    ax2[nt,1].scatter(df.lavg2,df.stdlog2, c='purple')
    #annotate sidf large figs
    ax2[nt,0].text(1.2,0.5,'Buffer: '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)

    #stdv vs mean

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.show()
fig1.savefig('../figures/dynamics_buffers_miliQ.png')
fig2.savefig('../figures/logDynamics_buffers_miliQ.png')



print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
