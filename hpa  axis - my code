# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 21:55:39 2026

@author: sathm
"""

# -*- coding: utf-8 -*-


import math as mt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
%matplotlib 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import find_peaks
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
    Ki1=Ki         #1*(1+kt);


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
    '''
    if inputF == "no":
        if t >= 98.25+d and t < 98.5+d or t >= 122.0+d and t < 122.25+d or t >=146.0+d and t < 146.25+d or t >= 170.0+d and t < 170.25+d or t >= 194.0+d and t < 194.25+d or t >= 218.0+d and t < 218.25+d or t >= 242.0+d and t < 242.25+d or t >= 266.0+d and t < 266.25+d:
            dex=dx*10;      #Comment this to run IC50 test
        #if t >= 98.0+d and t < 98.25+d or t >= 122.0+d and t < 122.25+d or t >=146.0+d and t < 146.25+d or t >= 170.0+d and t < 170.25+d or t >= 194.0+d and t < 194.25+d or t >= 218.0+d and t < 218.25+d or t >= 242.0+d and t < 242.25+d or t >= 266.0+d and t < 266.25+d: dex=dx*10;      #Comment this to run IC50 test

    '''
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

    # Cortisol - in cell
    dy4= (V_cort*StARp-k_cort*CORT);

    # Dexamethasone kinetics
    dy5= Dex-k_dex1*Dex1;
    dy6= V_dex2*Dex1-k_dex2*Dex2;

    # Peripheral cortisol - in plasma
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
    pq = 1; #do simulation with various values of pq 1,0.1
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

Ki=p_span3
n=p_span2
dex=p_span1
R=p_span4
Q=p_span5


ff=900
tspan=[0, ff];
kk= ff*10
t_eval = np.linspace(0, ff, kk) #from 0 to ff, with kk as step count



inputF = 'yes'

y0=[24.95694763498119,4.9123730894505435,0.9938828562096904,7.393415730911929,0.0,0.0,7.636200249054682,18.001529048278428,350.07408500237455,19.87971898349494,28.182728997762037,6.65145308687877e-05,41.97935885594562,50.96107608767422,0.2872674163673849,1.0053860612262895,3.5041839210619,1.072968514441532,7.636200249054682,1.0053860612262895,1.072968514441532,28.182728997762037,0.2872674163673849]
yd = solve_ivp(fun= Hpa_infla, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=y0 , rtol = 1e-10,atol = 1e-10)



#df_chk=pd.DataFrame(yd['y'][6]) #cortisol all timepoints values


#var=yd.y[6, 4000:] #var is the values of cortisol from 4000 onwards

#idx_max = var[np.argmax(var)] #highest value from that array

'''
#to find peaks, troughs and their indices for all 23 variables after reaching steady state (not just a single max or min)


t = yd.t                          # all 9000 time points in an array
timeframe = t[4000:]              # time points from 4000 onwards
y_frame = yd.y[:, 4000:]          # shape (23, len(timeframe))

# Peaks
peaks_indices = [find_peaks(y_frame[i])[0].tolist()   # shape[0] refers to the number of rows, [0] to only take the index
                 for i in range(y_frame.shape[0])]
peaks_values  = [y_frame[i, find_peaks(y_frame[i])[0]].tolist()
                 for i in range(y_frame.shape[0])]

# Troughs
troughs_indices = [find_peaks(-y_frame[i])[0].tolist()  # we use negative
                   for i in range(y_frame.shape[0])]
troughs_values  = [y_frame[i, find_peaks(-y_frame[i])[0]].tolist()
                   for i in range(y_frame.shape[0])]


peaks_times = [timeframe[x].tolist() for x in peaks_indices]
troughs_times = [timeframe[x].tolist() for x in troughs_indices] # true time from 0 to 900


# to check for cortisol

print(peaks_values[6])
print(peaks_times[6])


# to get a table for a certain value, just change the number

Peak_Values = peaks_values[6]
Peak_Times = peaks_times[6]
peak_df = pd.DataFrame({ "Peak Value": Peak_Values, "Peak Time": Peak_Times })
print(peak_df)

Trough_Values = troughs_values[6]
Trough_Times = troughs_times[6]
trough_df = pd.DataFrame({ "Trough Value": Trough_Values, "Trough Time": Trough_Times })
print(trough_df)

'''

'''
#to find amplitude, max peak, time of max peak of all 23 variables

results = []
for i in range(yd.y.shape[0]): #as it is a 2d array, we are iterating over the rows (each variable), hence we use shape[0]
    var = yd.y[i, :]
    amp = np.max(var) - np.min(var) #gives diff between max and min variable
    idx_max = np.argmax(var) #index of max peak
    idx_min = np.argmin(var) #index of min trough
    t_max = yd.t[idx_max] #what time point is in that index
    t_min = yd.t[idx_min] #what time point is in that index
    val_max = var[idx_max] #what variable is in that index
    val_min = var[idx_min] #what variable is in that index
    results.append((i+1, amp, val_max, t_max, val_min, t_min))

df = pd.DataFrame(results, columns=["Variable number", "Max-Min", "Max Peak", "Time of Max", "Min Trough", "Time of Min"])
print(df)

df["Amplitude"] = df["Max-Min"] / 2


'''

'''

#to create a plot with TNF, IL6, and CORTp

plt.plot(t_eval, yd['y'][6], label="CORTp", color="green") #dy7
plt.plot(t_eval, yd['y'][17], label="IL6", color="blue") #dy18
plt.plot(t_eval, yd['y'][15], label="TNF", color="red") #dy16

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.title("TNF, IL6, and Peripheral Cortisol")
plt.legend()

'''

'''
#to create two different subplots for pro-inflammatory and anti-inflammatory

f, (plt1, plt2) = plt.subplots(2) #creating diff subplots

line1 = plt1.plot(t_eval, yd['y'][17], label="IL6", color="blue") #pro inflammatory - IL6 and TNF
line1 = plt1.plot(t_eval, yd['y'][15], label="TNF", color="red")



line2 = plt2.plot(t_eval, yd['y'][16], label="IL10", color="blue") #anti inflammatory - IL10 and TGF
line2 = plt2.plot(t_eval, yd['y'][14], label="TGF", color="red")


#plt.plot(t_eval, yd['y'][17], label="IL6", color="blue") #dy18
#plt.plot(t_eval, yd['y'][15], label="TNF", color="red") #dy16
#plt.plot(t_eval, yd['y'][16], label="IL10", color="black")#dy17
#plt.plot(t_eval, yd['y'][14], label="TGF", color="green")#dy15



plt2.set_ylabel("Time")
plt2.set_xlabel("Concentration")
plt1.set_title("Pro-inflammatory vs Anti-inflammatory")
plt1.legend()
plt2.legend()

'''

'''
#to create diff subplots for TNF, TGF, IL6, and CORTp

f, (plt1, plt2, plt3, plt4, plt5) = plt.subplots(5) #creating diff subplots

line1 = plt1.plot(t_eval, yd['y'][17], label="IL6", color="blue") #dy18
#plot time points, at each timepoint
#saying we want every row in the first column to plot each timepoint
line2 = plt2.plot(t_eval, yd['y'][6]/10, label="CORTp", color="green") #dy7
line3 = plt3.plot(t_eval, yd['y'][14], label="TGF", color="red")#dy15
line4 = plt4.plot(t_eval, yd['y'][15], label="TNF", color="orange") #dy16
line5 = plt5.plot(t_eval, yd['y'][9], label="GRprotein", color="black") #dy10


plt4.set_xlabel('Time')
plt1.set_ylabel('IL6')
plt2.set_ylabel('CORTp')
plt3.set_ylabel('TGF')
plt4.set_ylabel('TNF')
plt5.set_ylabel('grprotein')

'''


plt.show()




