'''
All Hepes trials from Fig 1_a

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Goal is to have all Hepes data shown and then have analytical solutions for all that differ only in production (S_hooh)terms...with decay staying the same for all.

Using Morrisand Zinser 2013 Fig 1_a  data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_5_Hepes_Morris_2013/src/

to do: Get model working for Hepes treatments in data 
'''


import scipy 
from scipy.integrate import odeint
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   



#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 14})
plt.rcParams['legend.fontsize'] = 'small'
 


######################################################################

#  Importing the data as dfs for each trial

######################################################################

#/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/Project_5_Hepes_Morris_2013/data/inits

#import csv
df_all = pd.read_csv("../data/Hepes_assay_2013.csv") #use relative paths 
df_all.rename(columns = {'time(days)':'time'},inplace = True) 
df_all.fillna(0)
#use df clip? or just use 0H as control for base level of media 

#making log of data to look at error (how to do in loop? all combined whe I tried)
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




#calculating stats for trial compilation
df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4','rep5','rep6']]],axis=0)
df_all['sigma'] = np.nanstd(np.r_[[df_all[i] for i in ['rep1','rep2','rep3','rep4','rep5', 'rep6']]],axis=0)

df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3', 'log4','log5','log6']]],axis=0)
df_all['log_sigma'] = np.nanstd(np.r_[[df_all[i] for i in ['log1','log2','log3','log4','log5','log6']]],axis=0)


#split ROS treatments by number 
ROSs = df_all['Hepes_treatment'].unique()
ROSs = ROSs[~pd.isna(ROSs)] #clipping off nans 

nROSs = ROSs.shape[0]

colors = ('yellowgreen','green','c','b','purple','orange','r','k')
markers = ('.','P','s','1','^','*','p','+')


################
#graph for data viz 
#################

#fig creation and config 
fig1,ax1 = plt.subplots(figsize=[9,6])


plt.xlabel('Time (days)', fontsize = 20)
plt.ylabel('HOOH concentration (\u03BCM)', fontsize = 20)
plt.suptitle('HEPES Buffer HOOH Dynamics')
ax1.semilogy()
plt.tick_params(size = 14)


fig2,(ax2)= plt.subplots(nROSs,2, figsize = (13,11))
fig2.suptitle('Raw H dynamics of Hepes Buffer ')
fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.35)
ax2[0,0].set_title('Production Curve')
ax2[3,0].set_ylabel('HOOH concentration (\u03BCM)')
ax2[-1,0].set_xlabel('Time' )
ax2[0,0].semilogy()
ax2[0,1].set_title('STATS')
ax2[3,1].set_ylabel('STDV')
ax2[-1,1].set_xlabel('Mean' )

fig3,(ax3) = plt.subplots(nROSs,2,figsize = (13,11))
fig3.suptitle(' Log H dynamics of Hepes Buffer ')
fig3.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.35)
ax3[0,0].set_title('Production Curve')
ax3[3,0].set_ylabel('HOOH concentration (\u03BCM)')
ax3[-1,0].set_xlabel('Time' )
ax3[0,0].semilogy()
ax3[0,1].set_title('STATS')
ax3[3,1].set_ylabel('Log STDV')
ax3[-1,1].set_xlabel('Log Mean' )

fig4,(ax4) = plt.subplots(nROSs,2,figsize = (13,11))
fig4.suptitle(' Coorelations of Hepes ')
fig4.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.85, wspace=0.25, hspace=0.4)
ax4[0,0].set_title('Raw ')
ax4[0,1].set_title('Logged ')

#graphing ROS dfs each 
for (ros,ri) in zip(ROSs,range(nROSs)): # loop over ROS
    df = (df_all[df_all['Hepes_treatment']==ros]) # select one ROS treatment at a time 
    #ratios for graphing?
    raw_ratio =  (np.max(df.abundance/df.sigma))
    log_ratio =  (np.max(df.log_abundance/df.log_sigma))
    rats  = [raw_ratio,log_ratio]
    print(raw_ratio,log_ratio,rats)
    (rats) = np.rint(rats)
    best_ratio = (np.max(rats))
    rat_dif = (rats/best_ratio)
    stretch_rat = abs(rat_dif -1)
    print(stretch_rat)
    #plot
    ax1.plot(df.time,df.abundance,color = colors[ri], marker = markers[ri], markersize = 10, linestyle = ':', label='HEPES = '+str(ros)) # graph each 
    #raw data and stats for multi graph
    ax2[ri,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker = markers[ri],markersize= 10, label =('Mean'), c=colors[ri])
    ax2[ri,0].errorbar(df.time,df.avg1, yerr=df.stdv1,marker = markers[ri],markersize= 10, label =('avg1'), color = 'k' )
    ax2[ri,0].errorbar(df.time,df.avg2, yerr=df.stdv2,marker = markers[ri],markersize= 10, label =('avg2'), color = 'brown' )
    #log y axis and add legend to dynamics graph 
    ax2[ri,0].semilogy()

    l2 = ax2[ri,0].legend(loc = 'upper left')
    l2.draw_frame(False)
    #graph loggedstats
    ax2[ri,1].scatter(df.abundance,df.sigma, c=colors[ri])
    ax2[ri,1].scatter(df.avg1,df.stdv1, c='k')
    ax2[ri,1].scatter(df.avg2,df.stdv2, c='brown')
    
    #annotate side of large fig
    ax2[ri,0].text(1.2,0.5,('Hepes [] : '+ str(ros)),horizontalalignment='center', verticalalignment='center', transform=ax2[ri,1].transAxes)

###########
        #graph logged dynamics 
     
    ax3[ri,0].errorbar(df.time,df.log_abundance, yerr=df.sigma, marker = markers[ri],markersize= 10, label =('Log Mean'), c=colors[ri] )
    ax3[ri,0].errorbar(df.time,df.lavg1, yerr=df.stdlog1, marker = markers[ri],markersize= 10, label =('Log avg1'), color = 'k' )
    ax3[ri,0].errorbar(df.time,df.lavg2, yerr=df.stdlog2,marker = markers[ri],markersize= 10, label =('Log avg2'), color = 'brown' )
    #log y axis and add legend to dynamics graph 
    ax3[ri,0].semilogy()
    l3 = ax3[ri,0].legend(loc = 'upper left')
    l3.draw_frame(False)
    #annotate large graph 
    ax3[ri,0].text(1.2,0.5,'Hepes [] : '+ str(ros),horizontalalignment='center', verticalalignment='center', transform=ax3[ri,1].transAxes)
    
    #lgraph ogged stats
    ax3[ri,1].scatter(df.log_abundance,df.log_sigma, c=colors[ri])
    ax3[ri,1].scatter(df.lavg1,df.stdlog1, c='k')
    ax3[ri,1].scatter(df.lavg2,df.stdlog2, c='brown')


#create coorelations for fig 4
    raw_r,raw_p  = scipy.stats.pearsonr(df['abundance'],df['sigma'])
    log_r,log_p = scipy.stats.pearsonr(df['log_abundance'],df['log_sigma'])
    print(raw_r,log_r)
#Graph raw and logged coorelations for each ROS 


    #annotate coor graph 
    ax4[ri,0].hist(raw_r,color = 'b')
    ax4[ri,1].hist(log_r,color = 'r')
    ax4[ri,0].text(1.2,0.5,'Hepes [] : '+ str(ros),horizontalalignment='center', verticalalignment='center', transform=ax3[ri,1].transAxes)
 

    #working to set y and x lims
    '''
    #ax2.set_xlim(x.min()/ymin(), x.max()/ymax())
    #ax2[ri,1].set_xlim(0, int(np.max()))
    #plt.xlim(0, best_ratio)
    #plt.ylim(0, best_ratio)
    '''

#config large stats  graph labels 
# ylabels
for a in ax2[:,1]:
    a.set_ylabel('STDV')
    a.set_xlabel('Mean')
    #a.set_xlim(((stretch_rat[0])*df.abundance[0]),(stretch_rat[0])*df.abundance[-1])
    #a.set_ylim()

for a in ax3[:,1]:
    a.set_ylabel('Log STDV')
    a.set_xlabel('Log Mean ')
    #a.set_xlim()
    #a.set_ylim()



#makinng legend (must be after graphing to have label handles)
l1 = ax1.legend(loc = 'lower right', prop={"size":13}) 
l1.draw_frame(False)#print(df)
fig1.subplots_adjust(right=0.90, wspace = 0.25, hspace = 0.30) #shift white space for better fig view

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)


plt.show()
fig1.savefig('../figures/Hepes_all_data')
fig2.savefig('../figures/logDynamics_hepes.png')
fig3.savefig('../figures/rawdynamics_hepes.png')    
    
    

    #inits = pd.read_csv(("../data/inits/hepes"+str(ros)+".csv"))  













'''
####################################

 # Set Up ---OLD _______ model parts 

####################################
    
step = 0.05 #delta t
ndays = 26
times = np.linspace(0,ndays,int(ndays/step))

h_cons = np.array([0.15,0.6,1.5,1.5])
#h_convert = 0.6     #term used to convert input Hepes concentration to HOOH felt by cells
S_Hs = np.array([])

#making HOOH [ ] fromo HEPES [ ] given 
for r in ROSs: 
    S_HOOH = r * h_convert      #TAKING HEPES AT EACH TREATMENT AND FINDING HOOH BASED ON COMMON CONVERSION FACTOR
    S_Hs = np.append(S_Hs,S_HOOH)

for (h_con,ros) in zip(h_cons,ROSs): # loop over ROS
    S_HOOH = (ros * h_con)
    print(ros,S_HOOH)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    S_Hs = np.append(S_Hs,S_HOOH)

nS_Hs = S_Hs.shape[0]
deltas = np.array([0.06,0.04,0.03,0.01])

#model for use in ode int
def HsODEint(y,t,S_HOOH,delta):
    H = y[0]
    dHdt = S_HOOH-delta*H
    #print(dHdt,t,H)
    return dHdt 

liROS = ROSs.tolist()

#delta = 0.02

for (ros,S_HOOH,h0) in zip(ROSs,S_Hs,h0s): # loop over ROS
    count = (liROS.index(ros))
    print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,h0,times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')

for (ros,S_HOOH,delta) in zip(ROSs,S_Hs,deltas): # loop over ROS
    count = (liROS.index(ros))
    #print(ros,h0)   #the 1.0 and 0.1 start Hs are switched. Not sure where this occcurs. 
    solutions = odeint(HsODEint,r_[[0.7]],times,args=(S_HOOH,delta))  
    plt.plot(times,solutions[:,0],color = colors[count], marker = 'None')

plt.xlim([0, 26])
plt.ylim([0.15, 50])


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



plt.legend()
#plt.show()

fig, (ax2,ax3) = plt.subplots(1,2)
fig.suptitle('HEPES conversions',fontsize = 20)
plt.subplots_adjust(wspace = 0.3, top = 0.85)

#f2, (ax2, ax3) = plt.subplots(1,2, figsize=[8,6])
ax2.plot(ROSs,S_Hs, marker = 's', markersize = 8, color = 'brown')
ax2.set_title('HEPES vs modeled S_HOOH',fontsize = 10)
ax2.set_xlabel('HEPES added', fontsize = 10)
ax2.set_ylabel('Modeled HOOH Supply', fontsize = 10)

ax3.plot(ROSs,h_cons, marker = 'd', markersize = 8, color = 'k')
ax3.set_title('h_convert',fontsize = 10)
ax3.set_xlabel('HEPES added', fontsize = 10)
ax3.set_ylabel('h_convert needed for model', fontsize = 10)
#plt.tick_params(labelsize = 18) 
#plt.tick_params(size = 14)
#plt.legend(loc='lower right')

fig, (ax4) = plt.subplots()
fig.suptitle('Delta dynamics',fontsize = 20)
plt.subplots_adjust(wspace = 0.3, top = 0.85)


ax4.plot(ROSs,deltas, marker = 'd', markersize = 8, color = 'k')
#ax3.set_title('delta',fontsize = 10)
ax4.set_xlabel('HEPES added', fontsize = 10)
ax4.set_ylabel('delta needed for model', fontsize = 10)

'''
print('\n ~~~****~~~****~~~ \n')
print('Done with this H E P E S')
print('\n ~~~****~~~****~~~ \n')