import numpy as np
import MDAnalysis as MD
import os

# directories 
zero='/lustre/project/m2_trr146/asourpis/paper_ccn_h2o_0.25_0.75/75/NVTrun100ns'
uzero=MD.Universe(''+str(zero)+'/sim_100ns_0.0-V_nm.tpr',''+str(zero)+'/sim_100ns_0.0-V_nm.xtc')

#CHARGES
n1=uzero.select_atoms('resid 10 and name N1').total_charge()
c2=uzero.select_atoms('resid 10 and name C2').total_charge()
c1=uzero.select_atoms('resid 10 and name C1').total_charge()
h1=uzero.select_atoms('resid 10 and name H1').total_charge()
h2=uzero.select_atoms('resid 10 and name H2').total_charge()
h3=uzero.select_atoms('resid 10 and name H3').total_charge()

ow=uzero.select_atoms('resid 1500 and name OW').total_charge()
hw=uzero.select_atoms('resid 1500 and name HW1').total_charge()

#the calculation of the dot product mM for cental H2O collecting CCN and water after
def gwc(c,r,t):
    #central
    datac=np.loadtxt('/lustre/project/m2_trr146/asourpis/last_dance/kirkwood_pub/center_h2o/one_mol_h2o/r20pnts/better_stat_dt10/mol_1220/central_'+str(c)+'h2o_ccn_r'+str(r)+'_t'+str(t)+'.xvg',comments=['#','$','@'])[1::]
    testc=np.split(datac,len(datac)/3)
    ph2o=testc[0]*ow+testc[1]*hw+testc[2]*hw
    mi=ph2o
    #collection
    if os.path.isfile('/lustre/project/m2_trr146/asourpis/last_dance/kirkwood_pub/center_h2o/one_mol_h2o/r20pnts/better_stat_dt10/mol_1220/collect_'+str(c)+'h2o_ccn_r'+str(r)+'_t'+str(t)+'.xvg') == True:
        datacl=np.loadtxt('/lustre/project/m2_trr146/asourpis/last_dance/kirkwood_pub/center_h2o/one_mol_h2o/r20pnts/better_stat_dt10/mol_1220/collect_'+str(c)+'h2o_ccn_r'+str(r)+'_t'+str(t)+'.xvg',comments=['#','$','@'])[1::]
        fl=np.split(datacl,len(datacl)/18)
        dt1=[np.split(fl[i],6) for i in range(len(fl))]
        pccn=[dt1[j][0]*n1+dt1[j][1]*c2+dt1[j][2]*c1+dt1[j][3]*h1+dt1[j][4]*h2+dt1[j][5]*h3 for j in range(len(dt1))]
        M=sum(pccn)
    else:
        M=mi
    return np.dot(mi,M),np.dot(mi,mi)

def gww(c,r,t):
    #central
    datac=np.loadtxt('/lustre/project/m2_trr146/asourpis/last_dance/kirkwood_pub/center_h2o/one_mol_h2o/r20pnts/better_stat_dt10/mol_1220/central_'+str(c)+'h2o_h2o_r'+str(r)+'_t'+str(t)+'.xvg',comments=['#','$','@'])[1::]
    testc=np.split(datac,len(datac)/3)
    ph2o=testc[0]*ow+testc[1]*hw+testc[2]*hw
    mi=ph2o
    #collection
    if os.path.isfile('/lustre/project/m2_trr146/asourpis/last_dance/kirkwood_pub/center_h2o/one_mol_h2o/r20pnts/better_stat_dt10/mol_1220/collect_'+str(c)+'h2o_h2o_r'+str(r)+'_t'+str(t)+'.xvg') == True:
        datacl=np.loadtxt('/lustre/project/m2_trr146/asourpis/last_dance/kirkwood_pub/center_h2o/one_mol_h2o/r20pnts/better_stat_dt10/mol_1220/collect_'+str(c)+'h2o_h2o_r'+str(r)+'_t'+str(t)+'.xvg',comments=['#','$','@'])[1::]
        fl=np.split(datacl,len(datacl)/9)
        dt1=[np.split(fl[i],3) for i in range(len(fl))]
        ph2o=[dt1[j][0]*ow+dt1[j][1]*hw+dt1[j][2]*hw for j in range(len(dt1))]
        M=sum(ph2o)
    else:
        M=mi
    return np.dot(mi,M),np.dot(mi,mi)
    
c=1220
def Gofrccn(r,t):    
    mMccn=np.array(gwc(c,r,t)[0]/gwc(c,0.5,t)[1])
    return mMccn
    
c=1220
def Gofrh2o(r,t):    
    mMh2o=np.array(gww(c,r,t)[0]/gww(c,0.5,t)[1])
    return mMh2o    
    
tinit=1
tfnl=25000
dt=10
def kirkw(r):
    stath2o=[Gofrh2o(r,ts) for ts in range(tinit,tfnl,dt)]
    statccn=[Gofrccn(r,ts) for ts in range(tinit,tfnl,dt)]
    Nt=len(statccn)
    return np.average(stath2o),np.average(statccn),np.std(stath2o)/np.sqrt(Nt),np.std(statccn)/np.sqrt(Nt),np.std(statccn)+np.std(stath2o)/np.sqrt(Nt)


rls=[0.01, 0.26, 0.51, 0.76, 1.01, 1.26, 1.51, 1.76, 2.01, 2.26, 2.51, 2.76, 3.01, 3.26, 3.51, 3.76, 4.01, 4.26, 4.51, 4.76]
kirk=[kirkw(str(r)) for r in rls]
    
yh2o=[kirk[i][0] for i in range (len(kirk))]
yccn=[kirk[i][1] for i in range (len(kirk))]
yerrorh2o=[kirk[i][2] for i in range (len(kirk))]
yerrorccn=[kirk[i][3] for i in range (len(kirk))]
def norm():
    fnl=[]
    for i in range (len(kirk)):
        if yh2o[i]==yccn[i]:
            y=yh2o[i]
        else:
            y=yh2o[i]+yccn[i]
        fnl.append(y)
    return fnl
ytot=norm()
yerrortot=[kirk[i][4] for i in range (len(kirk))]

np.savetxt('gkh2o',yh2o)
np.savetxt('dgkh2o',yerrorh2o)
np.savetxt('gkccn',yccn)
np.savetxt('dgkccn',yerrorccn)
np.savetxt('gktot',ytot)
        
