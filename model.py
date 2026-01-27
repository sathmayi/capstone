# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 10:47:01 2022

@author: user
"""
import math as mt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
%matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random
import pandas as pd



def Hpa_infla(t,y):
    CRH=y[0];
    ACTH=y[1];
    StARp=y[2];
    CORT=y[3];
    Dex1 =y[4]; 
    Dex2 =y[5];

    CORTp =y[6];
    
    GR_mrna =y[7];
    GR_prot =y[8];
    GR_cyt =y[9];
    GR =y[10];

    LPS=y[11];
    Phg=y[12];
    Phg1= y[13];
    TGF=y[14];
    TNF=y[15];
    IL10=y[16];
    IL6=y[17];
    CORTp2 =y[18]
    TNFs = y[19]
    IL6s=y[20]
    GRs=y[21]
    TGFs=y[22]



    #Specify parameter perturbations here
    t1=0;
    te=5000;
    st=0            #*((t>t1) & (t<te));
    jt=r            #*((t>t1) & (t<te));
    kt=0            #*r*((t>t1) & (t<te));
    n=1             #*(1+jt);
    nj=n;           #1*(1+jt);
    Ki1=Ki;         #1*(1+kt);
    
    
    #Specify LPS dose here
    lp=p;
    It=lp
    #if t > 711 and t < 711.1:    It=10e7 *0.4        #10e7*0.4 #*10e7*((t>712) & (t<712.1));# for with cir3 10e7 *0.09 
    
    #uncomment this to disconnect HPA axis and inflammation to perform IC50 test
    #s0=0;
    
    
    #Specifiy dexamethasone dose here
    dx=1;
    dex=0;      #Comment this to run IC50 test
    #d=0;
    if inputF == "no":
        if t >= 98.25+d and t < 98.5+d or t >= 122.0+d and t < 122.25+d or t >=146.0+d and t < 146.25+d or t >= 170.0+d and t < 170.25+d or t >= 194.0+d and t < 194.25+d or t >= 218.0+d and t < 218.25+d or t >= 242.0+d and t < 242.25+d or t >= 266.0+d and t < 266.25+d:
            dex=dx*10;      #Comment this to run IC50 test 
        #if t >= 98.0+d and t < 98.25+d or t >= 122.0+d and t < 122.25+d or t >=146.0+d and t < 146.25+d or t >= 170.0+d and t < 170.25+d or t >= 194.0+d and t < 194.25+d or t >= 218.0+d and t < 218.25+d or t >= 242.0+d and t < 242.25+d or t >= 266.0+d and t < 266.25+d: dex=dx*10;      #Comment this to run IC50 test 
    Dex=dex           #*((t>696.00+d) & (t<696.25+d));
    
    
    # HPA-Axis
    K_strs=10;
    Ki_GR1=1.2        #*Ki ;
    n1=1              #*n;
    n2=2;
    k_crh=0.096;         
    V_acth1=1.4634;        #0.9756*1.5;
    k_acth=0.4312;         #1.15*0.1*2.5*1.5;          
    k_cort=0.99;           #0.44*0.9*2*1.25;     
    V_cort=5.67;           #0.0945*2*30;
    k_dex1=1.95;
    V_crh=0.0667           #*(1/R);
    V_acth2=11.2;
    Km_tnf3=40              #*R;
    V_starp2=15;
    k_starp=0.45;           #%0.3*1.5;
    V_starp1=0.0225;        #%0.015*1.5;
    Kil6=100;
    k_dex2=0.25;
    V_dex2=1.15;
    Ki_GR=50;
    V_GRmrna=3.4;
    km_GR1=25                   #*Ki;
    kon=0.00329;
    k_GRprot=0.0572 ;
    krt=0.63 ;
    kre=0.57;
    k_GRmrna=0.1124;
    V_GRprot=1.2;
    fr=0.49;
    tc=0.1401;
    ki_tnf=1e3;
    n3=2                             #*n;
    n4=2;
    fac=16.65;                  #0.5*33.3;
    
    #Circadian drive from SCN clock
    
    omega=2*mt.pi/24;
    
    cir=2*(1+ mt.cos(omega*(t)));
    
    #hh.append(cir)
    
    #Circadian cdrive from adrenal peripheral clock
    cir2=(4-cir)*(Ki_GR/(Ki_GR+GR));
    
    # Corticotrophic releasing hormone (CRH)
    dy1=(K_strs*(1+st)*(cir)*(Ki_GR1**n1/(Ki_GR1**n1+GR**n1))*(1+V_crh*TNF)-k_crh*CRH); #(K_strs*(1+st)*(cir)*(Ki_GR1**n1/(Ki_GR1**n1+GR**n1))*(1+V_crh*TNFs)-k_crh*CRH)
    
    # Adrenocorticotropic hormone (ACTH)
    dy2= (V_acth1*CRH*(Ki_GR1**n1/(Ki_GR1**n1+GR**n1))*(1+V_acth2*((TNF)**n2/(Km_tnf3**n2+ TNF**n2)))-k_acth*ACTH) ; #(V_acth1*CRH*(Ki_GR1**n1/(Ki_GR1**n1+GR**n1))*(1+V_acth2*((TNFs)**n2/(Km_tnf3**n2+ TNFs**n2)))-k_acth*ACTH)
    
    # StAR proetin
    dy3=s0*(V_starp1*(ACTH*cir2)*(1+V_starp2*(TNF**n2/(Km_tnf3**n2+TNF**n2))*(Kil6/(Kil6+IL6)))  -k_starp*StARp); #*cir2 #s0*(V_starp1*(ACTH*cir2)*(1+V_starp2*(TNFs**n2/(Km_tnf3**n2+TNFs**n2))*(Kil6/(Kil6+IL6s)))  -k_starp*StARp)
    
    # Cortisol
    dy4= (V_cort*StARp-k_cort*CORT);
    
    # Dexamethasone kinetics
    dy5= Dex-k_dex1*Dex1;
    dy6= V_dex2*Dex1-k_dex2*Dex2;
    
    # Peripheral cortisol
    dy7=(1/(tc))*(CORT+(fac*Dex2)-CORTp);
    
    # GR mRNA
    dy8 = V_GRmrna*(1 - (GR**n3/((km_GR1)**n3+GR**n3)))- k_GRmrna*GR_mrna; 
    
    # GR protein
    dy9 = V_GRprot*GR_mrna + fr*kre*GR- kon*(CORTp)*GR_prot - k_GRprot*GR_prot;
    
    # GR-cortisol complex in cytosol
    dy10= kon*(CORTp)*GR_prot- krt*GR_cyt*(ki_tnf**n4/(TNF**n4+ki_tnf**n4)); #kon*(CORTp)*GR_prot- krt*GR_cyt*(ki_tnf**n4/(TNFs**n4+ki_tnf**n4))
    
    # Nuclear GR
    dy11= krt*GR_cyt*(ki_tnf**n4/(TNF**n4+ki_tnf**n4)) -kre*GR; #krt*GR_cyt*(ki_tnf**n4/(TNFs**n4+ki_tnf**n4)) -kre*GR
    
    #cortisol conc in prepheral without dex
    dy19=(1/(tc))*(CORT-CORTp2);
    
    
    '''%%%%% Inflammatory pathway%%%%%%%%%#
    
    #%%%%%%%Validation parametset 2 for 0.4 ng/kg @ 2 PM%%%%%%%%%%%%%'''
    
    k_lps= 2.7e-5;         #1.35*1e-7*100*2;%0
    V_phg1= 4.9956*1e7*0.35 #4.9956*1e7*0.25 ;
    V_phg2=12.949;   #12.949
    Km_tnf1=1693.9509;
    Ki_tgf1=0.00721;
    Ki_IL10=7.384;        #147.68*0.05;
    k_phg=1.439;
    V_tgf1=0.15625*1e-8;
    k_tgf=0.0635;           #0.03177*2;
    
    V_tnf1=25.5194;  
    Km_phg1=412500;         #0.075*550*1e4;
    Ki_tgf2=0.143;          #0.15893*0.9;
    V_tnf2=106542;          # 3.0*3.5514*1e4;
    Km_tnf2=123.96;         #0.08*1.5495*1e3;
    k_tnf=1.25;
    Km_il61=80;
    
    Km_phg2=161012;          # 8.0506*1e7*0.2*1e-2;
    k_il10=1.6;              #200*0.008;
    V_il102=2.1938e3;         #43875*0.05;
    Km_tgf3=0.76;             # 0.38*2;
    V_il101=1.3374e3;          #2.67480*1e6*50.0e-5;
    Ki_il102=23.636;           #1.1818*20*1;
    V_il62=5.0e5;
    Km_il6tnf=339.164;        #4.8452*70;
    Km_phg3=11e6;              #110*10e4;
    V_il61= 0.55e5 
    k_il6=1.625;              #0.5*3.25;
    tp=1.5;
    tnf_b=1.25;
    il6_b=1.5;
    V_tgf2=0.5;
    Km_GR2=500                  #*Ki1;
    n5=1                         #*nj;
    n6=4;
    n7=2;
    n8=2;
    n9=6;
    n10=2;
    pq = 1; #do ssimulation with various values of pq 1,0.1
    cir3=4*(1+ mt.cos(omega*(t)+0.5))
    
    # LPS
    dy12= 1e-7*(1+It)-k_lps*LPS*Phg;
    
    # Phagocytes
    
    dy13= (V_phg1*(((1*cir3)+ (V_phg2* TNF/(Km_tnf1+TNF)))* (Ki_tgf1/(Ki_tgf1+ TGF))* (Ki_IL10/(Ki_IL10+IL10)))*LPS- k_phg*Phg);
    #dy13= (V_phg1*((1+ (V_phg2* TNF/(Km_tnf1+TNF)))* (Ki_tgf1/(Ki_tgf1+ TGF))* (Ki_IL10/(Ki_IL10+IL10)))*LPS- k_phg*Phg);
    
    dy14=(1.0/(tp))*(Phg-Phg1);
    
    # Transformig grwth factor (TGF)
    dy15=  (V_tgf1* Phg + (V_tgf2*GR**n5/(Km_GR2**n5+GR**n5))- k_tgf*TGF);
    
    # Tumour necrosis factor (TNF)
    dy16= (tnf_b+(Phg/(Km_phg1+Phg))*(V_tnf1+V_tnf2*(TNF/(Km_tnf2+TNF)))*(Ki_tgf2**n6/(Ki_tgf2**n6+TGF**n6))*(1-(IL6**n7/(Km_il61**n7+IL6**n7)))-k_tnf*(TNF))
    
    # Interleukin 10 (IL10)
    dy17=((V_il101*Phg1**n8/(Phg1**n8+Km_phg2**n8))+ (V_il102* TGF**n9/(Km_tgf3**n9+ TGF**n9)) - k_il10*IL10);
    
    # Interleukin 6 (IL6)
    dy18=(il6_b+(Phg1/(Km_phg3+Phg1))*(V_il61+V_il62*(TNF+IL6)**n10/(Km_il6tnf**n10+(TNF+IL6)**n10))*(Ki_il102/(Ki_il102+(IL10)))-k_il6*IL6)
    
    
    #TNFs
    dy20= (TNF - TNFs)*pq
    
    #IL6
    dy21= (IL6-IL6s)*pq
    
    #GR effect
    dy22= (GR - GRs)*pq
    #TGF slower
    dy23=(TGF-TGFs)*pq
    
    return [dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16,dy17,dy18,dy19,dy20,dy21,dy22,dy23]
    
#function for extracting average max and min from the range dynamics of particular bio-molecules
def max_avg(a,index):
       max1=[]
       n=0
       
       
       while n<len(a) :
           max1.append(max(a[n:n+240]))
           n+=241
       #print(max(a))
       #return [(sum(max1)/len(max1)),round(np.float64(t_eval[np.where(norm_yd['y'][index]== max(a))]),2)]
       return [(sum(max1)/len(max1)),round(np.float64(t_eval[np.where(yd['y'][index]== max(a))]),2)]
def min_avg(a,index):
       min1=[]
       n=0
       while n<len(a) :
           min1.append(min(a[n:n+240]))
           n+=241
       #return [(sum(min1)/len(min1)),round(np.float64(t_eval[np.where(norm_yd['y'][index]== min(a))]),2)]
       return [(sum(min1)/len(min1)),round(np.float64(t_eval[np.where(yd['y'][index]== min(a))]),2)]


# function for converting the time into AM-PM format
def time_conver(time,n):
    t1={}
    s1=0
    n1=n
    chk=0
    for i in time:
        frac, whole = mt.modf(i)
        if whole == n1:
            
            if frac >= 0.0 and frac < 0.25:
                s1=1
                t1.update({round(i, 2):str(s1)+" AM"})
            elif round(frac, 2) >= 0.25 and round(frac, 2) <= 0.50:
                t1.update({round(i, 2):str(s1)+":15 AM"})
            elif round(frac, 2) >= 0.50 and round(frac, 2) <= 0.75:
                t1.update({round(i, 2):str(s1)+":30 AM"})
            elif round(frac, 2) >= 0.75 and round(frac, 2) <= 0.80:
                t1.update({round(i, 2):str(s1)+":45 AM"})
            elif round(frac, 2) >= 0.80 and round(frac, 2) <= 0.99:
                    t1.update({round(i, 2):str(s1+1)+" AM"})
                    
                    
        elif i>n1 and i<(n1+12):
            if chk != whole:
                s1+=1
                chk=whole
            if round(frac, 2) >= 0.25 and round(frac, 2) <= 0.50 and s1 !=12:
                t1.update({round(i, 2):str(s1)+":15 AM"})
            elif round(frac, 2) >= 0.50 and round(frac, 2) <= 0.75 and s1 != 12:
                t1.update({round(i, 2):str(s1)+":30 AM"})
            elif round(frac, 2) >= 0.75 and round(frac, 2) <= 0.80 and s1!=12:
                t1.update({round(i, 2):str(s1)+":45 AM"})
            elif s1 == 12:
                if round(frac, 2) >= 0.25 and round(frac, 2) <= 0.50:
                    t1.update({round(i, 2):str(s1)+":15 PM"})
                elif round(frac, 2) >= 0.50 and round(frac, 2) <= 0.75:
                    t1.update({round(i, 2):str(s1)+":30 PM"})
                elif round(frac, 2) >= 0.75 and round(frac, 2) <= 0.80:
                    t1.update({round(i, 2):str(s1)+":45 PM"})
                elif round(frac, 2) >= 0.80 and round(frac, 2) <= 0.99:
                    #print(whole)
                    t1.update({round(i, 2):str(1)+" PM"})
                    
                else:
                    t1.update({round(i, 2):str(s1)+" PM"})
               
            elif frac>= 0.80 and frac<= 0.99:
                if s1==11:
                    t1.update({round(i, 2):str(s1+1)+" PM"})
                else:
                    t1.update({round(i, 2):str(s1+1)+" AM"})
            else:
                t1.update({round(i, 2):str(s1)+" AM"})   
     
        
        elif i>=(n1+12) and i<(n1+24):
            if chk != whole:
                s1+=1
                chk=whole
                #print(chk)
            
            if round(frac, 2) >= 0.25 and round(frac, 2) <= 0.50 and s1!=24:
                t1.update({round(i, 2):str(s1-12)+":15 PM"})
            elif round(frac, 2) >= 0.50 and round(frac, 2) <= 0.75 and s1!=24:
                t1.update({round(i, 2):str(s1-12)+":30 PM"})
            elif round(frac, 2) >= 0.75 and round(frac, 2) <= 0.80 and s1!=24:
                t1.update({round(i, 2):str(s1-12)+":45 PM"})
                
            elif s1 == 24:
                
                if round(frac, 2) >= 0.25 and round(frac, 2) <= 0.50:
                    t1.update({round(i, 2):str(s1-12)+":15 AM"})
                elif round(frac, 2) >= 0.50 and round(frac, 2) <= 0.75:
                    t1.update({round(i, 2):str(s1-12)+":30 AM"})
                elif round(frac, 2) >= 0.75 and round(frac, 2) <= 0.80:
                    t1.update({round(i, 2):str(s1-12)+":45 AM"})
                elif round(frac, 2) >= 0.80 and round(frac, 2) <= 0.99:
                    t1.update({round(i, 2):str(1)+" AM"})
                    if(round(frac, 2) > 0.90 and round(frac, 2) <= 0.99):
                        n1=n1+24
                        #print(s1)
                        s1=0
                        continue
                else:
                    t1.update({round(i, 2):str(s1-12)+" AM"})
                
                
            elif round(frac, 2) >= 0.80 and round(frac, 2) <= 0.99 and s1!=24:
                if s1==23:
                    t1.update({round(i, 2):str(s1-11)+" AM"})
                else:
                    t1.update({round(i, 2):str(s1-11)+" PM"})
            else:
                t1.update({round(i, 2):str(s1-12)+" PM"})    
            
               

    #print (n)
    return t1

##############################

p=0.0; #LPS injection
r=1.0; # Change this to vary GR sensitivity
p1=0;  # Change this to 1 to introduce LPS for IC50 test
s0=1;  # Change this to 0 to disconnect HPA with Inflmmation for IC50 test
#%u=1;  # 0.1% make this 0.1 for simulationg IC50 else 1

# Varying dexamethasone

p_span1=1;        #[0,logspace(-2,2.5,20)];%[0:0.01:0.2];%%[5 10 50 100];%[5 10 20 50];% 1.5 2 2.5 3 3.5 4];10;%
# Varying sensitivity
p_span2=1;        #[0.25:0.25:2];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];
# Varying inhibitory constant
p_span3=1;        #[0.1:0.1:1.5];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];

# Varying cytokine effect on HPA axis
p_span4=1;        #[0.1:0.1:1.5];%[0.25:0.25:2];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];

# Varying GR negative feedback
p_span5=1;       #[0.1:0.1:2];%;%[0.5:0.5:4];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];


# The data are extarcted from Grigolet et al., 2010, Clodie et al., 2008., Copeland et al., 2005,
# Lauw et al., 2005, Wegner et al., 2017

# Immune response for 0.4 ng/kg LPS%%%% 
# The data is extracted from Grigoleit et al.,2010.
T1=[0, 1,1.5 ,2,3,4,6];
il6x=[3.3,13.2,63.7,132,83.5,30.8,7.7];
tnfx=[9,72,89,51,16.5,12.4,10];
tnfs=[2.33,7.7,10.7,6.1,6.5,3.7,2.33];
il10x=[8.5,17.9,30.8,36.8,34.6,15.4,10.7];

T2=[0,1.5,3,6];
cortx=[10,11,19,6];
Acthx=[8.0, 11.0,25.0,5.0];
T3=[0,1.5,2,3,4,5,6];
cortx1=[8,10,16,19,18,15,12];

'''
% %%%Immune response for 2 ng/kg LPS%%%% 
% T2=[0, 1,1.5 ,2,3,4,6];
% il6x1=[6.86, 96,528,1010,761,322,20.6];
% tnfx1=[3.32,267,445,302,109,26.1,20];
% il10x1=[11.9,42.9,81,176,260,114,33.3];
% cortx2=[10.1,9.08,15,17,20,19,13];
% Acthx2=[32,31,40,75,74,35,22];

for count3=1:1:length(p_span3)
     Ki=p_span3(count3);
     
for count2=1:1:length(p_span2)
n=p_span2(count2);

 for count1=1:1:length(p_span1)
 dex=p_span1(count1);
 
 for count4=1:1:length(p_span4)
     R=p_span4(count4);
     
for count5=1:1:length(p_span5)
     Q=p_span5(count5);
'''

Ki=p_span3
n=p_span2
dex=p_span1
R=p_span4
Q=p_span5 


ff=900
tspan=[0, ff];
kk= ff*10
t_eval = np.linspace(0, ff, kk)

t_conver= time_conver(t_eval,98.0)

# Initial conditions
#y0=[17.4896364862443,4.82352867112865,0.794555311165795,10.8943173854454,0,0,11.1479571299593,20.1160322148203,342.786973832252,24.8798181706611,31.3302385380763,0.000129677801623061,5897.54196590606,25,0.446309929849742,1.00160367257638,22.8093965060946,1.09564141359877];
#y0=[24.95694763498119,4.9123730894505435,0.9938828562096904,7.393415730911929,0.0,0.0,7.636200249054682,18.001529048278428,350.07408500237455,19.87971898349494,28.182728997762037,6.65145308687877e-05,41.97935885594562,50.96107608767422,0.2872674163673849,1.0053860612262895,3.5041839210619,1.072968514441532]

#gg1=[46.7199320899054,25.389872802971016,3.630239339118392,19.830336583555372,0.0,0.0,19.654680820529816,25.113101554794866,374.3426144370158,29.410114259658794,20.52784653524815,6.661147881474037e-05,82.1296395378398,79.49978017569399,0.19137911992054096,1.0355290873813494,0.3242202136173553,1.1591583045576057]


#yd = solve_ivp(fun= Hpa_infla, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=y0)

max_cort_dex=[]
max_cort=[]
max_tnf=[]
max_il6=[]
min_cort=[]
min_cort_dex=[]
min_tnf=[]
min_il6=[]
cort_deff=[]
aftr_dex=[]
PhaseDiff_cort_dex=[] 
PhaseDiff_cort=[] 
PhaseDiff_tnf=[]
PhaseDiff_il6=[]
n_change=0
   


F=1
norm_yd=[]

inputF = input("Do u want to run the normal dynamics also Answer in yes or no > ")
if inputF == "yes":
    y0=[24.95694763498119,4.9123730894505435,0.9938828562096904,7.393415730911929,0.0,0.0,7.636200249054682,18.001529048278428,350.07408500237455,19.87971898349494,28.182728997762037,6.65145308687877e-05,41.97935885594562,50.96107608767422,0.2872674163673849,1.0053860612262895,3.5041839210619,1.072968514441532,7.636200249054682,1.0053860612262895,1.072968514441532,28.182728997762037,0.2872674163673849]
    yd = solve_ivp(fun= Hpa_infla, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=y0 , rtol = 1e-10,atol = 1e-10)
    norm_yd=yd
    F=0
    #inputF ="no"

inputF ="no"   

#print(y0[14])

'''
for i in np.arange(0,26,):
    d=i
    y0=[24.95694763498119,4.9123730894505435,0.9938828562096904,7.393415730911929,0.0,0.0,7.636200249054682,18.001529048278428,350.07408500237455,19.87971898349494,28.182728997762037,6.65145308687877e-05,41.97935885594562,50.96107608767422,0.2872674163673849,1.0053860612262895,3.5041839210619,1.072968514441532,7.636200249054682,1.0053860612262895,1.072968514441532,28.182728997762037,0.2872674163673849]
    yd = solve_ivp(fun= Hpa_infla, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=y0 , rtol = 1e-10,atol = 1e-10)
    #plt.plot(t_eval,norm_yd['y'][15],'r')
    #plt.plot(t_eval,yd['y'][15],'b')
    #plt.plot(t_eval,yd['y'][16], 'r') 
    #plt.show()
    max_cort.append(max_avg(yd['y'][18][1530:2661+n_change],18))
    max_cort_dex.append(max_avg(yd['y'][6][1530:2661+n_change],6))
    max_tnf.append(max_avg(yd['y'][15][1530:2661+n_change],15))
    max_il6.append(max_avg(yd['y'][17][1530:2661+n_change],17))
    min_cort.append(min_avg(yd['y'][18][1530:2661+n_change],18))
    min_cort_dex.append(min_avg(yd['y'][6][1530:2661+n_change],6))
    min_tnf.append(min_avg(yd['y'][15][1530:2661+n_change],15))
    min_il6.append(min_avg(yd['y'][17][1530:2661+n_change],17))
    PhaseDiff_cort_dex.append(float(t_eval[np.where(yd['y'][6]== max(yd['y'][6][2660:2900]))]-t_eval[np.where(norm_yd['y'][6] == max(norm_yd['y'][6][2660:2900]))]))
    PhaseDiff_cort.append(float(t_eval[np.where(yd['y'][18]== max(yd['y'][18][2660:2900]))]-t_eval[np.where(norm_yd['y'][18] == max(norm_yd['y'][18][2660:2900]))]))
    PhaseDiff_tnf.append(float(t_eval[np.where(yd['y'][15]== max(yd['y'][15][2660:2900]))]-t_eval[np.where(norm_yd['y'][15] == max(norm_yd['y'][15][2660:2900]))]))
    PhaseDiff_il6.append(float(t_eval[np.where(yd['y'][17]== max(yd['y'][17][2660:2900]))]-t_eval[np.where(norm_yd['y'][17] == max(norm_yd['y'][17][2660:2900]))]))
    n_change=10
  

#print(len(norm_yd['y'][16][980:4820]))
print(t_conver)

#tconver2=str(t_conver.values())


a=max_avg(norm_yd['y'][16][1530:1991],16)[0]
b=max_avg(norm_yd['y'][16][1530:1991],16)[1]
print('il10',a,t_conver[b]) # taking out the max peak time
del a
del b
a=min_avg(norm_yd['y'][16][1530:1991],16)[0]
b=min_avg(norm_yd['y'][16][1530:1991],16)[1]
print('il10',a,t_conver[b]) # taking out the max peak time
del a
del b



plt.plot(t_eval,norm_yd['y'][8], 'r')

plt.plot(t_eval[1530:1991],norm_yd['y'][18][1530:1991], 'g')

plt.xlabel("Time", fontweight='bold')
plt.ylabel("Conc",fontweight='bold')
plt.legend()
#plt.xticks(tconver2, fontweight='bold')
plt.show()





lp1,lp2,lp3,lp4,lp5,lp6=norm_yd['y'][18][7100:7190],norm_yd['y'][1][7100:7190],norm_yd['y'][14][7100:7190],norm_yd['y'][15][7100:7190],norm_yd['y'][16][7100:7190],norm_yd['y'][17][7100:7190]


plt.subplot(3,2,1)
plt.plot(t_eval[7100:7190],norm_yd['y'][18][7100:7190], 'bo', label= "Cort with cir3",markersize=3.5) 
plt.plot(t_eval[7100:7190],lp1, 'm-.',linewidth=2, label= "Cort w/o cir3") 
plt.legend()
plt.subplot(3,2,2)
plt.plot(t_eval[7100:7190],norm_yd['y'][1][7100:7190], 'bo', label= "ACTH with cir3",markersize=3.5)
plt.plot(t_eval[7100:7190],lp2, 'm-.',linewidth=2, label= "ACTH w/o cir3")  
plt.legend()
plt.subplot(3,2,3)
plt.plot(t_eval[7100:7190],norm_yd['y'][14][7100:7190], 'bo', label= "TGF with cir3",markersize=3.5)
plt.plot(t_eval[7100:7190],lp3, 'm-.',linewidth=2, label= "TGF w/o cir3")  
plt.legend()
plt.subplot(3,2,4)
plt.plot(t_eval[7100:7190],norm_yd['y'][15][7100:7190], 'bo', label= "TNF with cir3",markersize=3.5)
plt.plot(t_eval[7100:7190],lp4, 'm-.',linewidth=2, label= "TNF w/o cir3")  
plt.legend()
plt.subplot(3,2,5)
plt.plot(t_eval[7100:7190],norm_yd['y'][16][7100:7190], 'bo', label= "IL10 with cir3",markersize=3.5) 
plt.plot(t_eval[7100:7190],lp5, 'm-.',linewidth=2, label= "IL10 w/o cir3") 
plt.legend()
plt.subplot(3,2,6)
plt.plot(t_eval[7100:7190],norm_yd['y'][17][7100:7190], 'bo', label= "IL6 with cir3",markersize=3.5) 
plt.plot(t_eval[7100:7190],lp6, 'm-.',linewidth=2, label= "IL6 w/ocir3") 
plt.legend()






fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(t_eval,norm_yd['y'][6], 'r') 
ax.grid()
plt.plot(t_eval,hh,'r')
plt.plot(t_eval,hh3,'b')
plt.xlabel("Time", fontweight='bold')
plt.ylabel("Conc",fontweight='bold')
plt.show()



gh=pd.DataFrame({'alpha':, ' '})

gh.to_csv


ampli_cort=(np.subtract([s[0] for s in max_cort],[s[0] for s in min_cort]))
ampli_cort_dex=(np.subtract([s[0] for s in max_cort_dex],[s[0] for s in min_cort_dex]))
ampli_tnf=(np.subtract([s[0] for s in max_tnf],[s[0] for s in min_tnf]))
ampli_il6=(np.subtract([s[0] for s in max_il6],[s[0] for s in min_il6]))

norm_cort=(max(norm_yd['y'][18][1530:2661])-min(norm_yd['y'][18][1530:2661]))
norm_tnf=(max(norm_yd['y'][15][1530:2661])-min(norm_yd['y'][15][1530:2661]))
norm_il6=(max(norm_yd['y'][17][1530:2661])-min(norm_yd['y'][17][1530:2661]))
norm_cort_dex=(max(norm_yd['y'][6][1530:2661])-min(norm_yd['y'][6][1530:2661]))

graph= pd.DataFrame({'Cort_max':[s[0] for s in max_cort],'Cort_min':[s[0] for s  in min_cort],"PhaseDiff-Cort":PhaseDiff_cort, "ampitude":ampli_cort})
graph1=pd.DataFrame({'TNF_max':[s[0] for s in max_tnf], 'TNF_min':[s[0] for s in min_tnf], "PhaseDiff-TNF":PhaseDiff_tnf, "ampitude":ampli_tnf})
graph2=pd.DataFrame({'IL6_max':[s[0] for s in max_il6], 'IL6_min':[s[0] for s in min_il6], "PhaseDiff-IL6":PhaseDiff_il6, "ampitude":ampli_il6})
graph3= pd.DataFrame({'Cort-Dex_max':[s[0] for s in max_cort_dex],'Cort-Dex_min':[s[0] for s  in min_cort_dex],"PhaseDiff-CortDex":PhaseDiff_cort_dex, "ampitude":ampli_cort_dex})



graph.to_csv("D:\samaya\work_related\graph_CortDayna.csv")
graph1.to_csv("D:\samaya\work_related\graph_TNFDayna.csv")
graph2.to_csv("D:\samaya\work_related\graph_IL6Dayna.csv")
graph3.to_csv("D:\samaya\work_related\graph_CortDexDayna.csv")


graph4= pd.DataFrame({'Time':t_eval[980:1220],'TNF_comp-dynamics':norm_yd['y'][15][980:1220]})

graph4.to_csv("D:\samaya\work_related\graph_normTNF.csv")
del graph4
del graph1
del graph2
del graph3
del graph

  

plt.xticks(xnumbers, fontweight='bold')
plt.yticks(fontweight='bold')
ax.grid()
#plt.axis([1,25,1.11,1.16]) # [xstart, xend, ystart, yend]
plt.xlabel("Time(AM-PM)", fontweight='bold')
plt.ylabel("Cort_max_Conc",fontweight='bold')
plt.legend()
plt.show()




#ax.plot(xnumbers,[s[0] for s in max_cort], 'g',linewidth=2, linestyle ='-')  
for i in range(0,9):
  ax.text(xnumbers[i], max_cort[i][0],t_conver[max_cort[i][1]], size=9,fontstretch = 'extra-condensed', fontstyle ='italic',fontweight="bold",verticalalignment ="top" ,color ="b" )
for i in range(9,25):
  ax.text(xnumbers[i], max_cort[i][0],t_conver[max_cort[i][1]], size=9,fontstretch = 'extra-condensed', fontstyle ='italic',fontweight="bold",verticalalignment ="bottom" ,color='b')

plt.xticks(xnumbers, fontweight='bold')
plt.yticks(fontweight='bold')
ax.grid()
#plt.axis([1,25,1.11,1.16]) # [xstart, xend, ystart, yend]
plt.xlabel("Time(AM-PM)", fontweight='bold')
plt.ylabel("Cort_max_Conc",fontweight='bold')
plt.legend()
plt.show()




xnumbers = range(1, 26, 1)
plt.plot(xnumbers,PhaseDiff_cort, 'r')  
plt.xticks(xnumbers)
plt.grid()
plt.xlabel("Time(AM-PM)")
plt.ylabel("Phase diff il6 ")
plt.legend()
plt.show()


hh=[]
hh1=[]
hh3=[]
p=0
for i in t_eval :
    omega=2*mt.pi/24;
    hh.append(2*(1+ mt.cos(omega*(i))))
for k in t_eval:
    hh3.append(8*(1+ mt.cos(omega*(k)+0.5)))
    
for j in norm_yd['y'][10]:
    hh1.append((4-hh3[p])*(50/(50+j)))
    p+=1

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(t_eval,norm_yd['y'][6], 'r') 
ax.grid()
#ax.plot(t_eval,norm_yd['y'][15], 'r') 
#plt.plot(t_eval,norm_yd['y'][18], 'r', label= "") 
#plt.plot(t_eval,norm_yd['y'][17], 'b') 
#plt.plot(t_eval,hh1,'r', label= "Cir2")
plt.plot(t_eval,hh3,'b', label= "Cir3")
plt.plot(t_eval,hh,'g', label= "Cir1")
plt.xlabel("Time", fontweight='bold')
plt.ylabel("impulse",fontweight='bold')
plt.legend()
plt.show()

'''


'''
  
d=12
y0=[24.95694763498119,4.9123730894505435,0.9938828562096904,7.393415730911929,0.0,0.0,7.636200249054682,18.001529048278428,350.07408500237455,19.87971898349494,28.182728997762037,6.65145308687877e-05,41.97935885594562,50.96107608767422,0.2872674163673849,1.0053860612262895,3.5041839210619,1.072968514441532]
yd = solve_ivp(fun= Hpa_infla, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=y0 , rtol = 1e-10,atol = 1.e-10)

#max_cort.append(max(yd['y'][6][980+n_change:2661+n_change]))
plt.plot(t_eval,yd['y'][15], 'r') 
#plt.plot(t_eval,yd['y'][5], 'r')  
#plt.plot(t_eval,yd['y'][6], 'r')  
#plt.plot(t_eval,yd['y'][5], 'r') 
#plt.plot(t_eval,yd['y'][16], 'r') 
plt.show()
plt.xlabel("Angle in Radians")
plt.ylabel("Magnitude")
plt.title("Plot of some trigonometric functions")
plt.yticks(ynumbers)
plt.legend()


#for random color selection
  r1=random.randint(0,255)
  g1=random.randint(0,255)
  b1=random.randint(0,255)
  plt.plot(t_eval,yd['y'][17], color= (r1/255,g1/255,b1/255))
  plt.show()


#chk1=yd['y'][6]
plt.plot(t_eval,yd['y'][6], 'r') 
plt.plot(t_eval,chk1) 
plt.show()



y
gg=[]
for i in range(len(yd['y'])):
    gg.append(yd['y'][i][99931])
    
    
    
xnumbers = range(1, 26, 1)
plt.plot(xnumbers,PhaseDiff_tnf, 'r')  
plt.xticks(xnumbers)
plt.grid()
#plt.axis([1,25,1.11,1.16]) # [xstart, xend, ystart, yend]
plt.xlabel("Time(AM-PM)")
plt.ylabel("phase diff cort ")
plt.legend()
plt.show()



data={'Cort_max':[s[0] for s in max_cort],'Cort_min':[s[0] for s  in min_cort],"PhaseDiff-Cort":PhaseDiff_cort}
data1={'TNF_max':[s[0] for s in max_tnf], 'TNF_min':[s[0] for s in min_tnf], "PhaseDiff-TNF":PhaseDiff_tnf}
data2={'IL6_max':[s[0] for s in max_il6], 'IL6_min':[s[0] for s in min_il6], "PhaseDiff-IL6":PhaseDiff_il6}



graph= pd.DataFrame(data)
graph1=pd.DataFrame(data1)
graph2=pd.DataFrame(data2)

graph.to_csv("graph.csv")
graph1.to_csv("D:\samaya\work_related\graph1.csv")
graph2.to_csv("D:\samaya\work_related\graph2.csv")



a1,a2,a3,a4,a5,a6= norm_yd['y'][18],norm_yd['y'][8],norm_yd['y'][14],norm_yd['y'][15],norm_yd['y'][16],norm_yd['y'][17]  # without

#ax.plot(t_eval,norm_yd['y'][15], 'r') 

fig, ax = plt.subplots(3,2, figsize=[11, 9])




#plt.subplot(3,2,1)
ax[0,0].plot(t_eval[980:1700],norm_yd['y'][18][980:1700], 'r', label= "cort") 
ax[0,0].axvspan(98,103, color='gray')
ax[0,0].axvspan(116,126, color='gray')
ax[0,0].axvspan(140,151, color='gray')
ax[0,0].axvspan(164,170, color='gray')

#plt.plot(t_eval,a1, 'g')
ax[0,0].legend()

#plt.subplot(3,2,2)
ax[0,1].plot(t_eval[980:1700],norm_yd['y'][8][980:1700], 'b', label= "GR_prt") 
ax[0,1].axvspan(98,103, color='gray')
ax[0,1].axvspan(116,126, color='gray')
ax[0,1].axvspan(140,151, color='gray')
ax[0,1].axvspan(164,170, color='gray')

#plt.plot(t_eval,a2, 'g')
ax[0,1].legend()

#plt.subplot(3,2,3)
ax[1,0].plot(t_eval[980:1700],norm_yd['y'][14][980:1700], 'g', label= "TGF") 
ax[1,0].axvspan(98,103, color='gray')
ax[1,0].axvspan(116,126, color='gray')
ax[1,0].axvspan(140,151, color='gray')
ax[1,0].axvspan(164,170, color='gray')

#plt.plot(t_eval,a3, 'g')
ax[1,0].legend()

#plt.subplot(3,2,4)
ax[1,1].plot(t_eval[980:1700],norm_yd['y'][15][980:1700], 'y', label= "TNF") 
ax[1,1].axvspan(98,103, color='gray')
ax[1,1].axvspan(116,126, color='gray')
ax[1,1].axvspan(140,151, color='gray')
ax[1,1].axvspan(164,170, color='gray')

#plt.plot(t_eval,a4, 'g')
ax[1,1].legend()

#plt.subplot(3,2,5)
ax[2,0].plot(t_eval[980:1700],norm_yd['y'][16][980:1700], 'c', label= "IL10") 
ax[2,0].axvspan(98,103, color='gray')
ax[2,0].axvspan(116,126, color='gray')
ax[2,0].axvspan(140,151, color='gray')
ax[2,0].axvspan(164,170, color='gray')

#plt.plot(t_eval,a5, 'g')
ax[2,0].legend()

#ax[2,1].subplot(3,2,6)
ax[2,1].plot(t_eval[980:1700],norm_yd['y'][17][980:1700], 'm', label= "IL6") 
ax[2,1].axvspan(98,103, color='gray')
ax[2,1].axvspan(116,126, color='gray')
ax[2,1].axvspan(140,151, color='gray')
ax[2,1].axvspan(164,170, color='gray')
#plt.plot(t_eval,a6, 'g')
ax[2,1].legend()


'''
'''
# Adding a plot in the figure which will encapsulate all the subplots with axis showing only
fig.add_subplot(1, 1, 1, frame_on=False)

# Hiding the axis ticks and tick labels of the bigger plot
plt.tick_params(labelcolor="none", bottom=False, left=False)

# Adding the x-axis and y-axis labels for the bigger plot
plt.xlabel("Time",fontsize=10, fontweight='bold')
plt.ylabel("Concentration",fontsize=10, fontweight='bold')

plt.show()

'''




    
    
    
    
