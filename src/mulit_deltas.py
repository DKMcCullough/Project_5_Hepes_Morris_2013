'''
All Hepes trials from Fig 1_a

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Goal is to have all Hepes data shown and then have analytical solutions for all that differ only in production (S_hooh)terms...with decay staying the same for all.

Using Morrisand Zinser 2013 Fig 1_a  data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_5_Hepes_Morris_2013/scripts/
'''


from scipy import *
from scipy.integrate import odeint
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




######################################################################

#  Importing the data as dfs for each trial

######################################################################

#import csv
df_all = pd.read_csv("../data/fig1a_reformat.csv") #use relative paths 
df_all.rename(columns = {'Time (days)':'times','treatment (milliMolar)':'treatment', 'HOOH (micromolar)':'HOOH'}, inplace = True)
#df_all[['rep1', 'rep2','rep3']] = df_all[['rep1', 'rep2','rep3']].fillna(value=0)
#df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first
#df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'
#neeed to cut off extra NAN columns 


ROSs = df_all['treatment'].unique()
nROSs = ROSs.shape[0]
#ROSs.sort()     #Change here to get ambient ROS to be 61 or 90 for species and 0 is detoxed. 


f1,ax1 = plt.subplots(figsize=[8,6])
colors = ('green','b','orange','r')

h0s = np.array([])

for (ros,ri) in zip(ROSs,range(nROSs)): # loop over ROS
    tdf = (df_all[df_all['treatment']==ros]) # select one ROS treatment at a time 
    hoohs = tdf['HOOH']
    h0 = hoohs.iloc[0]
    h0s = np.append(h0,h0s) #fing start H for each trial for easier modeling later
    ax1.plot(tdf['times'],tdf['HOOH'], marker='s',color = colors[ri], markersize = 10, linestyle = ':', label='HEPES = '+str(ros)) # graph each 



plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel('HOOH concentration (\u03BCM)', fontsize = 20)
plt.yscale('log')
plt.tick_params(labelsize = 18) 
plt.tick_params(size = 14)
plt.legend(loc='upper left')
#plt.show()



####################################

 # Set Up model parts 

####################################
    
step = 0.05 #delta t
ndays = 26
times = np.linspace(0,ndays,int(ndays/step))

#h_cons = np.array([0.15,0.6,1.1,1.1])
h_convert = 0.6     #term used to convert input Hepes concentration to HOOH felt by cells
S_Hs = np.array([])

#making HOOH [ ] fromo HEPES [ ] given 
for r in ROSs: 
    S_HOOH = r * h_convert      #TAKING HEPES AT EACH TREATMENT AND FINDING HOOH BASED ON COMMON CONVERSION FACTOR
    S_Hs = np.append(S_Hs,S_HOOH)
'''
for (h_con,ros) in zip(h_cons,ROSs): # loop over ROS
    S_HOOH = ros * h_con
    print(ros,S_HOOH)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    S_Hs = np.append(S_Hs,S_HOOH)
'''
nS_Hs = S_Hs.shape[0]
deltas = np.array([0.4,0.04,0.0,0.0])

#model for use in ode int
def HsODEint(y,t,S_HOOH,delta):
    H = y[0]
    dHdt = S_HOOH-delta*H
    #print(dHdt,t,H)
    return dHdt 

liROS = ROSs.tolist()

#delta = 0.02
'''
for (ros,S_HOOH,h0) in zip(ROSs,S_Hs,h0s): # loop over ROS
    count = (liROS.index(ros))
    print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,h0,times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')
'''
for (ros,S_HOOH,delta) in zip(ROSs,S_Hs,deltas): # loop over ROS
    count = (liROS.index(ros))
    #print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,r_[[0.7]],times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')

plt.xlim([0, 26])
plt.ylim([0.15, 50])

'''
####################################

#analytical solutions 

####################################


S_HOOH = S_Hs

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 

plt.plot(times,Hs,c='y',marker='.',label='Analytical Solution')



delta = 0.17
for (ros,S_HOOH) in zip(ROSs,S_Hs): # loop over ROS
    count = (liROS.index(ros))
    #print(count)
    solutions = Hs = f(times,S_HOOH,delta)
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')

'''

plt.legend()
plt.show()

fig, (ax3) = plt.subplots()
fig.suptitle('Delta dynamics',fontsize = 20)
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax3.plot(ROSs,deltas, marker = 'd', markersize = 8, color = 'k')
ax3.set_title('delta',fontsize = 10)
ax3.set_xlabel('HEPES added', fontsize = 10)
ax3.set_ylabel('delta needed for model', fontsize = 10)
#plt.tick_params(labelsize = 18) 
#plt.tick_params(size = 14)
#plt.legend(loc='lower right')





print('Done')