'''
All Hepes trials from Fig 1_a

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Goal is to have all Hepes data shown and then have analytical solutions for all that differ only in production (S_hooh)terms...with decay staying the same for all.

Using Morrisand Zinser 2013 Fig 1_a  data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_5_Hepes_Morris_2013/scripts/
'''


from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




######################################################################

#  Importing the data as dfs for each trial

######################################################################
'''
#####   10 Hepes
Hepes_10_df = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx', header = 2)   #importing data while ignoring 2 header rows...from first sheet of excel file
Hepes_10_df.columns = ['Time (days)', 'HOOH Concentration (\u03BCM)']     #renaming columns of dataframe so they are labeled appropriately 

######  1 Hepes
Hepes_1_df = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx','1 mM Hepes', header = 2)   #importing data for 1.0 trial 
Hepes_1_df.columns = ['Time (days)', 'HOOH Concentration (\u03BCM)']     

 
###### 0.1 Hepes

Hepes_01_df = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx','0.1mM Hepes', header = 2)   #importing data for 0.1 trial
Hepes_01_df.columns = ['Time (days)', 'HOOH Concentration (\u03BCM)']    

#######   0 Hepes

Hepes_0_df = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx','0mM Hepes', header = 2)   #importing data for 0 trial 
Hepes_0_df.columns = ['Time (days)', 'HOOH Concentration (\u03BCM)']  

'''
#import csv
df_all = pd.read_csv("../data/fig1a_reformat.csv") #use relative paths 
df_all.rename(columns = {'treatment (milliMolar)':'treatment', 'HOOH (micromolar)':'HOOH'}, inplace = True)
#df_all[['rep1', 'rep2','rep3']] = df_all[['rep1', 'rep2','rep3']].fillna(value=0)
#df_all = df_all.dropna(axis = 1)     #taking NaN columns off the end of df but had to fill rep 1 and 2 Nans first
df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'
#neeed to cut off extra NAN columns 


ROSs = df_all['treatment'].unique()
#ROSs.sort()     #Change here to get ambient ROS to be 61 or 90 for species and 0 is detoxed. 


f1,ax1 = plt.subplots()
for r in ROSs:
    if (df_all['treatment']==r):
        ax1[].plot(df_all['times'],df_all['HOOH'],marker='o',label='Hepes [ ]'+ str(r))

'''


plt.scatter(Hepes_10_df[['Time (days)']],Hepes_10_df[['HOOH Concentration (\u03BCM)']], marker = 'x', s=50, c = 'r', label = 'HOOH from 10\u03BCM Hepes')
plt.scatter(Hepes_1_df[['Time (days)']],Hepes_1_df[['HOOH Concentration (\u03BCM)']], marker = 'x', s=50, c = 'b', label = 'HOOH from 1.0\u03BCM Hepes')
plt.scatter(Hepes_01_df[['Time (days)']],Hepes_01_df[['HOOH Concentration (\u03BCM)']], marker = 'x', s=50, c = 'g', label = 'HOOH from 0.1\u03BCM Hepes')
plt.scatter(Hepes_0_df[['Time (days)']],Hepes_0_df[['HOOH Concentration (\u03BCM)']], marker = 'x', s=50, c = 'y', label = 'HOOH from 0\u03BCM Hepes')


plt.xlabel('Time (day)', fontsize = 16)
plt.ylabel('HOOH concentration (\u03BCM)', fontsize = 16)
plt.yscale('log')
plt.tick_params(labelsize = 12) 
#plt.show()



####################################

 # Units for all analytical solutions

####################################

delta = 1.6    
step = 0.05 #delta t
ndays = 26
times = np.linspace(0,ndays,int(ndays/step))




####################################

#analytical solution for 10

####################################

#initial values and creating time array

S_HOOH = 15

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


plt.plot(times,Hs,c='r',marker='.',label='10\u03BCM Hepes Analytical Solution')




####################################

#analytical solution for 1

####################################

#initial values and creating time array

S_HOOH = 7.8

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


plt.plot(times,Hs,c='b',marker='.',label='1\u03BCM Hepes Analytical Solution')




####################################

#analytical solution for 0.1

####################################

#initial values and creating time array

S_HOOH = 3.4

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


plt.plot(times,Hs,c='g',marker='.',label='0.1\u03BCM Hepes Analytical Solution')




####################################

#analytical solution for 00

####################################

#initial values and creating time array

S_HOOH = 1.0

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


plt.plot(times,Hs,c='y',marker='.',label='0\u03BCM Hepes Analytical Solution')






plt.legend()
plt.show()
'''