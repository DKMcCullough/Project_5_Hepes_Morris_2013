'''
Hepes_production_of_hooh.py

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Using Morrisand Zinser 2013 Fig 1_a  data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_5_Hepes_Morris_2013/scripts/
'''


from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




###############################.

#  Importing the data frame for 10mM HOOH

###############################

 
Hepes_1_df = pd.read_excel('../data/data_fig1a_Morris_and_Zinser_2013.xlsx','1 mM Hepes', header = 2)   #importing data while ignoring 2 header rows
Hepes_1_df.columns = ['Time (days)', 'HOOH Concentration (\u03BCM)']     #renaming columns of dataframe so they are labeled appropriately 



#print(Hepes_1_df)                             #printing df to make sure I have all the data and labels correct
#print(Hepes_1_df[['Time (days)']])                  #making sure I can call the columns by thier new names
#print(Hepes_1_df[['HOOH Concentration (\u03BCM)']])


plt.scatter(Hepes_1_df[['Time (days)']],Hepes_1_df[['HOOH Concentration (\u03BCM)']], marker = 'x', s=50, c = 'b', label = 'HOOH from 1\u03BCM Hepes')
plt.xlabel('Time (day $^{-1}$)', fontsize=16)
plt.ylabel('HOOH concentration (\u03BCM)',fontsize=16)
plt.yscale('log')
plt.tick_params(labelsize=12)
#plt.show()



####################################

#analytical solution

####################################

#initial values and creating time array

delta = .7
S_HOOH = 4.5
step = 0.05 #delta t
ndays = 26
times = np.linspace(0,ndays,int(ndays/step))

def f(t, S_HOOH, delta):
        H = (S_HOOH/delta)*(1-e**(-delta*t))
        return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


plt.plot(times,Hs,c='b',marker='.',label='1\u03BCM Hepes Analytical Solution')


plt.legend()
plt.show()

