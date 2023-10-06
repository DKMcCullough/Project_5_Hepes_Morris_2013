'''
viz_dark_buffers_2013.py

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


df_all = pd.read_csv('../data/dark_buffer_2013MorrisSIfig1a.csv')

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#had to float then int for some reason so that we could calc log after 
df_all.astype({"rep1":"float","rep2":"float","rep3":"float","rep4":"float","rep5":"float","rep6":"float"})
df_all.astype({"rep1":"int","rep2":"int","rep3":"int","rep4":"int","rep5":"int","rep6":"int"})


#df_all.mask(df_all < 0, 0)
#df_all['rep6'].clip(lower=0)



#creating log and stats for data 
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])
df_all['log4'] = np.log(df_all['rep4'])
df_all['log5'] = np.log(df_all['rep5'])
df_all['log6'] = np.log(df_all['rep6'])



df_all['avg1'] = np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['avg2'] = np.nanmean(np.r_[[df_all[i] for i in ['rep4', 'rep5', 'rep6']]],axis=0)
df_all['stdv1'] = np.std(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['stdv2'] = np.std(np.r_[[df_all[i] for i in ['rep4', 'rep5', 'rep6']]],axis=0)

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



treats  = df_all['Buffer'].unique()
ntreats = treats.shape[0]



##############################

#    graphing the data 

##############################

fig1,(ax1)= plt.subplots(ntreats,2, figsize = (12,10))
fig1.suptitle('Raw H dynamics of Dark Buffer Production ')
fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax1[0,0].set_title('Production Curve')
ax1[1,0].set_ylabel('HOOH concentration (\u03BCM)')
ax1[3,0].set_xlabel('Time' )
ax1[0,0].semilogy()
ax1[0,1].set_title('STATS')
ax1[1,1].set_ylabel('STDV')
ax1[3,1].set_xlabel('Mean' )

fig2,(ax2) = plt.subplots(ntreats,2,figsize = (12,10))
fig2.suptitle(' Log H dynamics of Dark Buffer Production ')
fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax2[0,0].set_title('Production Curve')
ax2[1,0].set_ylabel('HOOH concentration (\u03BCM)')
ax2[2,0].set_xlabel('Time' )
ax2[0,0].semilogy()
ax2[0,1].set_title('STATS')
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
    l1 = ax1[nt,0].legend(loc = 'upper left')
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
    l2 = ax2[nt,0].legend(loc = 'upper left')
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
fig1.savefig('../figures/dynamics_dark_buffers.png')
fig2.savefig('../figures/logDynamics_dark_buffers.png')







print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
