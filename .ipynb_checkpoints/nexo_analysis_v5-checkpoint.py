#!/usr/bin/env python
# coding: utf-8

# In[1]:


import import_ipynb
from nexoheader import *
import uproot3


import sys


# In[2]:


def BlankNexoDF_Sens():
    
    return pd.DataFrame(columns= ['energy', 'energy_mctruth_allLXe', 'energy_mctruth_inTPC',
       'r_max_simple', 'lower_z', 'upper_z', 'standoff', 'standoff_r',
       'standoff_z', 'passed_z_thresh', 'passed_xy_thresh', 'r_max_3d',
       'n_x_ch_abovenoise', 'n_y_ch_abovenoise', 'max_r_mctruth',
       'lower_z_mctruth', 'upper_z_mctruth', 'standoff_mctruth',
       'standoff_r_mctruth', 'standoff_z_mctruth', 'max_r_active_mctruth',
       'lower_z_active_mctruth', 'upper_z_active_mctruth',
       'standoff_active_mctruth', 'standoff_r_active_mctruth',
       'standoff_z_active_mctruth', 'm_nQ', 'm_nOPCal', 'm_nOPCollected',
       'm_DNNvalue', 'NESTBugFound', 'ratio', 'fGenX', 'fGenY', 'fGenZ',
       'Branch_ratio', 'Seed', 'Nevts'])


# In[3]:


def TotalSimCount(DF):

    try:
        return float(DF[["Seed","Nevts"]].groupby("Seed").mean().sum())
        
    except:
        print(-1)
        
def Gauss(x, A, B,C):
    y = A*np.exp( (-1*(x-B)**2)/2*C**2)
    return y


# In[4]:


def GetNexoDF(tgt,version="-v3"):
    sensfile="/home/wrkshp/nexo/slurm_data/sens/%s%s-sens_1-200.root"%(tgt,version)



    File=uproot.open(sensfile)
    keys=File.keys()

    F_keys=File[keys[0]]
    F_dict={}
    for i,j in enumerate(F_keys.keys()):
        F_dict[j]=F_keys[j].array()

    DF=pd.DataFrame.from_dict(F_dict)    
    return DF


# In[5]:


def GetNexosingleroot(DF):
    
    info=DF.iloc[0]
    
    sensfile="/home/wrkshp/nexo/slurm_data/sens/%s-%s-seed%s.root"%(info["Target"],info["Version"],info["Seed"])
    
        
    

    File=uproot.open(sensfile)
    keys=File.keys()

    F_keys=File[keys[0]]
    F_dict={}
    for i,j in enumerate(F_keys.keys()):
        F_dict[j]=F_keys[j].array()

    DF=pd.DataFrame.from_dict(F_dict)    
    DF["Target"]=info["Target"]
    return DF


# In[ ]:





# In[6]:



standoff=100

sscuts_array=["standoff>%i and m_DNNvalue > 0.850 "%(standoff),
    " and n_x_ch_abovenoise > 0 and n_y_ch_abovenoise > 0",
    " and passed_xy_thresh == 1 and passed_z_thresh== 1 " ,
    " and (m_nOPCal <1.077*m_nQ+313) ",
    " and (m_nOPCal >0.597*m_nQ-216) "]


#cuts_array=["standoff>-200 and  standoff<200"]


sscuts= "".join(sscuts_array)


TL_peak_cut_array = sscuts_array + ["and abs(energy - 2614) < 10"]
TL_peakcuts= "".join(TL_peak_cut_array)


Co_peak_cut_array = sscuts_array + ["and abs(energy - 1332.5) < 10"]
Co_peakcuts= "".join(Co_peak_cut_array)

Cs_peak_cut_array = sscuts_array + ["and abs(energy - 661.657) < 10"]
Cs_peakcuts= "".join(Cs_peak_cut_array)


# In[7]:


NDB=NexoDB()


# In[196]:


AA=NDB.iloc[0]


# In[197]:


NDB.query("Target == 'Tl208' and Version == 'v4' and TfileError == 0")


# In[87]:



    


# In[ ]:





# In[8]:


class SourceSim:
    def __init__(self, target,version="v3"):
        self.Target=target
        self.Version=version
        
        
    def ReadDB(self,DB=""):
        
        if type(DB) != type(pd.DataFrame):
            DB=NexoDB()
            
        
        self.TotalSim = DB.query("Target == '%s' and Version == '%s'"%(self.Target, self.Version))["NumEvents"].sum()
        self.SeedCount= len(np.array(DB.query("Target == '%s' and Version == '%s'"%(self.Target, self.Version))["Seed"]))
        self.MissingSeedCount =len(DB.query("Target == '%s' and Version == '%s' and TfileError == 'G4 does not exsits'"%(self.Target, self.Version)))
        self.DB =DB.query("Target == '%s' and Version == '%s'"%(self.Target, self.Version))
        self.SeedList=DB.query("Target == '%s' and Version == '%s' and RootCheck == True"%(self.Target, self.Version))["Seed"]
    def LoadRoot(self):
        
        DFs=[]
        for i,j in enumerate(self.SeedList):
            
            
            DFs.append(GetNexosingleroot(self.DB.query("Seed == '%s'"%(j))))
            #print(j, type(GetNexosingleroot(self.DB.query("Seed == '%s'"%(j)))))
            

        if len(DFs) > 1:
            DF=pd.concat(DFs)
        else:
            DF=BlankNexoDF_Sens()
            
            
        
        return DF
            
            
            
            
            
    


# In[ ]:





# In[199]:


Sim_Th228=SourceSim("Tl208")


# In[12]:


Tl208v4 = SourceSim("Tl208","v4")
Tl208v4.ReadDB()
Tl208v4.SeedList , Tl208v4.SeedCount


# In[13]:


Tl208v4.DB


# In[14]:



Sim_Th228.ReadDB()
Sim_Th228.SeedCount, Sim_Th228.TotalSim, Sim_Th228.MissingSeedCount
Th228_RF=Sim_Th228.LoadRoot()
Th228_RF


# In[29]:


Sim_Th228.SeedList


# In[19]:


chain=["Th228","Ra224","Pb212","Bi212","Tl208"]
version="v3"
if 1==1:
#def Get DecayChain(chain,version="v3"):

    Chainsources=[]
    Chainsource_list=[]
    for tgt in chain:
        Source = SourceSim(tgt,version)
        Source.ReadDB()
        SourceDF=Source.LoadRoot()
        print(len(SourceDF["energy"]))
        Chainsources.append(SourceDF)
        Chainsource_list.append(Source)
    ChainDF=pd.concat(Chainsources)
        
        




# In[205]:


Tl208v4.Version


# In[192]:


ChainDF.to_pickle("Th288Chain_1M.pkl")


# In[162]:


ChainDF.keys()


# In[163]:


len(ChainDF["energy"]),ChainDF["Branch_ratio"].sum(),len(ChainDF.query("Target == 'Tl208'"))


# In[26]:



fig, ax = plt.subplots(figsize=(8.5, 5))
fig, bx = plt.subplots(figsize=(8.5, 5))
fig, cx = plt.subplots(figsize=(8.5, 5))

nbins=3000
peakcuts=TL_peakcuts

A=[]
B=[]
C=[]
deposited_events=0
dep_nobr=0
GEvents_sum=0
GEvents_nobr=0


A.append(
    ax.hist(
            x=ChainDF["energy"]/1000
            ,bins=nbins
            #,linewidth=0.5
            ,range=[0,3]
            ,alpha=0.55
            #,histtype='step'
            #,color="darkgrey"
            ,label="%s"%("Th228 Chain"),
            
       # weights=ChainDF["Branch_ratio"]
            )
        )


for i,tgt in enumerate(chain):
    
    
    pltDF=ChainDF.query("Target == '%s'"%(tgt))
    if len(pltDF["energy"])==0:
        continue
        
    print("Target :",tgt, ", Deposited Events :", len(pltDF["energy"]), ", \t Weight with BR:",
          round(pltDF["Branch_ratio"].sum(),1), ",  \t Total Sim: ",Chainsource_list[i].TotalSim)
    deposited_events=deposited_events+round(pltDF["Branch_ratio"].sum(),1)
    dep_nobr=len(pltDF["energy"])+dep_nobr
    
    SSDF=pltDF.query(sscuts)
    PeakDF=pltDF.query(peakcuts)
    
    GEvents=SSDF["Branch_ratio"].sum()
    GEventsnbr=len(SSDF["energy"])
    
    GEvents_sum=GEvents_sum+GEvents
    GEvents_nobr=GEvents_nobr+GEventsnbr
    
    print("Alpha ",len(pltDF["energy"])/Chainsource_list[i].TotalSim )
    print("\t\tGEvents :", GEvents, ", Gamma ", GEvents/Chainsource_list[i].TotalSim,
         "GEvents nobr ", GEventsnbr, " Gammanbr ", GEventsnbr/Chainsource_list[i].TotalSim )
    
    A.append(
        ax.hist(
                x=pltDF["energy"]/1000
                ,bins=nbins
                ,range=[0,3]
                ,alpha=0.75
                ,histtype='step'
                #,   color=Tlcolor,
                ,label="%s"%(tgt),
        #        weights=pltDF["Branch_ratio"]
                )
            )
    B.append(
                bx.hist(
                x=SSDF["energy"]/1000
                ,bins=nbins
                ,range=[0,3]
                ,alpha=0.75
                ,histtype='step'
                #,   color=Tlcolor,
                ,label="%s"%(tgt),
         #       weights=SSDF["Branch_ratio"]
                )
    )
    
    
    C.append(
                cx.hist(
                x=(PeakDF["energy"]/1000)+(0.001*i)
                ,bins=nbins
                ,range=[0,3]
                ,alpha=0.75
                ,histtype='step'
                #,   color=Tlcolor,
                ,label="%s"%(tgt),
                weights=PeakDF["Branch_ratio"]
                )
    )
            
    
    
    

if 1==0:
    A.append(
        ax.hist(
            x=ChainDF["energy"]/1000
            ,bins=nbins
            ,linewidth=0.5
            ,range=[0,3]
            ,alpha=0.55
            ,histtype='step'
            ,color="darkgrey"
            ,label="%s"%("Th228 Chain"),
            
        weights=ChainDF["Branch_ratio"]
            )
        )    
    
#print("\n\nTh228 Chain ", " Deposited Events", deposited_events, " alpha :", deposited_events/Chainsource_list[0].TotalSim)

#print("Th228 Chain ", " Deposited Events no BR", dep_nobr, " alpha no BR :", dep_nobr/Chainsource_list[0].TotalSim)

#print("GEvents", GEvents_sum)
#print(" Gamma : ",GEvents_sum/Chainsource_list[0].TotalSim )
#print(" Gammanbr : ",GEvents_nobr/Chainsource_list[0].TotalSim )

#Gamma=GEvents_sum/Chainsource_list[0].TotalSim
#alpha=deposited_events/Chainsource_list[0].TotalSim
    




ax.legend(ncol=2, loc="upper center")
ax.set_yscale("log")
ax.set_ylabel("BR Weight \nEvents")
ax.set_xlabel("Energy(MeV)")
ax.grid(True)
ax.set_xlim(-0.2,3.15)

bx.legend(ncol=2, loc="upper center")
bx.set_yscale("log")
bx.set_ylabel("BR Weight \nEvents")
bx.set_xlabel("Energy(MeV)")
bx.grid(True)
bx.set_xlim(-0.2,3.15)

cx.legend(ncol=1, loc="upper right")
cx.set_ylabel("BR Weight \nEvents")
cx.set_xlabel("Energy(MeV)")
cx.set_xlim(2.4,2.9)
cx.grid(True)


# In[180]:


alpha=alpha
gamma=Gamma 
GoodEventReq=100

print(alpha, gamma)
t= 0.001
A=np.arange(1,50000,1)
  
f2 = np.exp(-alpha*A*t)
f3 =(1+alpha*A*t)*np.exp(-alpha*A*t)


rate_f3=gamma*A*f3
time_f3=GoodEventReq/rate_f3

rate_f2=gamma*A*f2
time_f2=GoodEventReq/rate_f2


timeto=min(GoodEventReq/(gamma*A*f3)/(3600))
act=(1+np.sqrt(5))/(2*alpha*t)

timeto_f2=min(GoodEventReq/(gamma*A*f2)/(3600))
act_f2=1/(alpha*t)

print( round(act_f2,3),round(timeto_f2,3))

fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(A/1000,y=rate_f2,label="Rate Th228 Chain ")

ax.set_ylabel("Rate(Hz)")
ax.set_xlabel("Activity(kBq)")


from scipy import interpolate
rate_est=interpolate.interp1d(A, rate_f2, kind='linear')

ax.grid("True")
ax.legend()
#ax.set_xscale("log")
print("Rate @ 4 kBq :", np.round(  rate_est(4000) ,3) , "Hz" )

print("Max Rate @ %.1f kBq -> %.3f Hz "%(act_f2, np.round(  rate_est(act_f2) ,3)))


# In[5]:


Tl208=GetNexoDF("Tl208")
Bi212=GetNexoDF("Bi212")


# In[ ]:




Tl208_totalsim=TotalSimCount(Tl208)
Bi212_totalsim=TotalSimCount(Bi212)


Tl208_alpha=len(Tl208.query("energy > 0")["energy"])/Tl208_totalsim
Tl208_gamma=len(Tl208.query(TL_peakcuts))/Tl208_totalsim
print("Tl208 :: alpha ",Tl208_alpha," Gamma ", Tl208_gamma)

Bi212_alpha=len(Bi212.query("energy > 0")["energy"])/Bi212_totalsim
Bi212_gamma=len(Bi212.query(TL_peakcuts))/Bi212_totalsim
print("Bi212 :: alpha ",Bi212_alpha," Gamma ", Bi212_gamma)


# In[ ]:





# In[ ]:


tgt="Cs137"
cut=Cs_peakcuts
version="-v3"





    
DF_totalsim=TotalSimCount(DF)

if DF_totalsim == -1:
    print("Issue with total sim calc")
    exit
    
    
DF_alpha=len(DF.query("energy > 0")["energy"])/DF_totalsim
DF_gamma=len(DF.query(cut))/Tl_totalsim


# In[ ]:


DF


# In[35]:



nbins=500

fig, ax = plt.subplots(figsize=(6, 4),dpi=300)
plt.hist(Tl208["energy"]/1000
         ,nbins
         ,range=[0,3]
         , alpha=0.75,
                    histtype='step',
                      #   color=Tlcolor,
                 label="All Deposits")

plt.hist(Tl208.query(sscuts)["energy"]/1000
         ,nbins
         ,range=[0,3]
         , alpha=0.75,
                    histtype='step',
                      #   color=Tlcolor,
                 label="SS")


plt.xlim(0,3)
plt.grid(True)
plt.legend(loc ="upper right")
plt.yscale("log")
plt.ylabel("Deposits")
plt.xlabel("Energy(MeV)")
#plt.title("Cs137 Deposits")


# In[ ]:




