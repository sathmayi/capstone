import math as mt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
%matplotlib
import matplotlib.pyplot as plt
import random
import pandas as pd



PP1=[]
def odes2(t,x):

    #print (t)
    per= x[0]
    cry= x[1]
    rev= x[2]
    ror= x[3]
    bmal= x[4]
    Prot_per= x[5]
    Prot_cry= x[6]
    Prot_rev= x[7]
    Prot_ror= x[8]
    Prot_bmal= x[9]
    PC= x[10]
    CB= x[11]
    nampt= x[12]
    Prot_nampt= x[13]
    NAD= x[14]
    dbp= x[15]

    ##parameters for mRNA and protein degradatation rates
    dm_bmal = 0.827333442085 #Bmal mRNA degradation rate
    dm_cry = 0.319848706181 #Cry mRNA degradation rate
    dm_dbp = 0.384486532062 #Dbp mRNA degradation rate
    dm_nampt = 0.814311309051 #Nampt mRNA degradation rate
    dm_per = 0.323114598647 #Per mRNA degradation rate
    dm_rev = 4.0479072846 #Rev-Erb mRNA degradation rate
    dm_ror = 0.246575760727 #Ror mRNA degradation rate
    dp_bmal = 0.186413466797 #BMAL protein degradation rate
    dp_cry = 0.599026119971 #CRY protein degradation rate
    dp_nampt = 49.8841023982 #NAMPT protein degradation rate
    dp_per = 10.9446515392 #PER protein degradation rate
    dp_rev = 0.281564063057 #REV-ERB protein degradation rate
    dp_ror = 0.0340112196281 #ROR degradation rate
    d_cb = 0.197714012552 #CLOCK-BMAL complex degradation rate
    d_pc = 0.609290138329 #PER-CRY complex degradation rate

    ##parameters for Complexation kinetic rates
    Kass_cb = 0.0162923358655 # nmol−1· l · h−1 CLOCK-BMAL association rate
    Kass_pc = 12.302005485 # nmol−1· l · h−1 PER-CRY association rate
    Kdiss_cb = 0.00613502224231 # h−1 CLOCK-BMAL dissociation rate
    Kdiss_pc = 0.0365867175408 # h−1 PER-CRY dissociation rate

    ##parmeters for Maximal transcription rates
    Vmax_bmal = 0.0916862770193 #Bmal maximal transcription rate
    Vmax_cry = 0.702216664807 #Cry maximal transcription rate
    Vmax_dbp = 0.0802009609453 #Dbp maximal transcription rate
    Vmax_nampt = 3.49035201578 #Nampt maximal transcription rate
    Vmax_per = 0.730201742662 #Per maximal transcription rate
    Vmax_rev = 1.12297601784 #Rev-Erb maximal transcription rate
    Vmax_ror = 6.9843472736 #Ror maximal transcription rate

    ##parameters for Activation ratios
    fold_bmal = 15.9258093373 #activation ratio of Bmal by ROR
    fold_cry = 1.1604489571 #activation ratio of Cry by CLOCK-BMAL
    fold_dbp = 400.0 #activation ratio of Dbp by CLOCK-BMAL
    fold_nampt = 1.57880681573 #activation ratio of Nampt by CLOCK-BMAL
    fold_per = 12.977351123 #activation ratio of Per by CLOCK-BMAL
    fold_rev = 73.2342431701 #activation ratio of Rev-Erb by CLOCK-BMAL
    fold_ror = 335.923333883 #activation ratio of Ror by CLOCK-BMAL

    ##parameters for
    Ka_bmal_ror = 0.00560038941378 #Regulation threshold of Bmal by ROR
    Ka_cry_cb = 1.0089387144 #Regulation threshold of Cry by CLOCK-BMAL
    Ka_dbp_cb = 0.308745016237 #Regulation threshold of Dbp by CLOCK-BMAL
    Ka_nampt_cb = 3.54586790835 #Regulation threshold of Nampt by CLOCK-BMAL
    Ka_per_cb = 2.03485134763 #Regulation threshold of Per by CLOCK-BMAL
    Ka_rev_cb = 0.260846828116 #Regulation threshold of Rev-Erb by CLOCK-BMAL
    Ka_ror_cb = 0.266407416327  #Regulation threshold of Ror by CLOCK-BMAL
    Ki_bmal_rev0 = 0.0108449480001 #Regulation threshold of Bmal by REV-ERB
    Ki_cry_rev0 = 0.248955507809 #Regulation threshold of Cry by REV-ERB
    Ki_cry_pc = 0.00338463577329 #Regulation threshold of Cry by PER-CRY
    Ki_dbp_pc = 2.23913672671 #Regulation threshold of Dbp by PER-CRY
    Ki_nampt_pc = 0.0137106537972 #Regulation threshold of Nampt by PER-CRY
    Ki_per_pc = 0.273493946059 #Regulation threshold of Per by PER-CRY
    Ki_rev_pc = 28.5885406354 #Regulation threshold of Rev-Erb by PER-CRY
    Ki_ror_pc = 0.0072858432208 #Regulation threshold of Ror by PER-CRY

    ##parameters for Hill coefficients
    hill11 = 4.32985205032 #Hill coeff., regulation of Bmal by REV-ERB
    hill10 = 1.83992599578 #Hill coeff., regulation of Bmal by ROR
    hill3 = 9.1109447538 #Hill coeff., regulation of Cry by CLOCK-BMAL
    hill4 = 2.43715119318 #Hill coeff., regulation of Cry by PER-CRY
    hill5 = 4.20952050286 #Hill coeff., regulation of Cry by REV-ERB
    hill14 = 7.32066818222 #Hill coeff., regulation of Dbp by CLOCK-BMAL
    hill15 = 10.4312927466 #Hill coeff., regulation of Dbp by PER-CRY
    hill12 = 1.91470474775 #Hill coeff., regulation of Nampt by CLOCK-BMAL
    hill13 = 1.34080593157 #Hill coeff., regulation of Nampt by PER-CRY
    hill1 = 8.52414053707 #Hill coeff., regulation of Per by CLOCK-BMAL
    hill2 = 8.53897990872 #Hill coeff., regulation of Per by PER-CRY
    hill6 = 9.83701536835 #Hill coeff., regulation of Rev-Erb by CLOCKBMAL
    hill7 = 3.31257899336 #Hill coeff., regulation of Rev-Erb by PER-CRY
    hill8 = 9.36456505724 #Hill coeff., regulation of Ror by CLOCK-BMAL
    hill9  = 1.84102439743 #Hill coeff., regulation of Ror by PER-CRY

    ##parameters for Translation rates
    Kp_bmal = 0.628507384997 #Bmal translation rate
    Kp_cry = 3.76912711677 #Cry translation rate
    Kp_nampt = 58.9647983852 #Nampt translation rate
    Kp_per = 13.2909782781 #Per translation rate
    Kp_rev = 0.0513221194724 #Rev-Erb translation rate
    Kp_ror = 0.0412765888526 #Ror translation rate

    ##parametrs for Protein stability modulation constants
    m_cry_ampk = 0.07940386211 #modulation of CRY stability by AMPK
    m_nampt_ampk = 0.6352243791 #modulation of NAMPT stability by AMPK
    m_per_ampk = 0.005243953731 #modulation of PER stability by AMPK
    m_per_sirt = 0.005452322816 #modulation of PER stability by SIRT

    ##parameters for NAD kinetics, Sirt1 and PGC1a activity
    Vsirt = 0.915854846427 #N/A Maximum SIRT1 activity
    Ksirt = 0.75 #nmol · l−1 Value of [NAD] at which SIRT1 activity is half of maximum
    d_nad = 378.009130749 #h−1 Rate of transformation of NAD into NAM
    Knad = 0.321746039086 #nmol · l−1 Value of [NAD] at which transformationinto NAM is at half of maximum rate
    NAD_basal = 0.9116166306 #nmol · l−1 Value of [NAD] below which transformation of NAD into NAM is inactivated
    Vnad = 296.3933869 #molecule per hour per NAMPT protein Maximum regeneration rate of NAD
    NAD_tot = 4.166901679 #nmol · l−1 Total concentration of NAD and NAM
    Knam = 2.76496 #nmol · l−1 Value of NAM at which NAD salvage rate is half of maximum
    Vpg = 24.06372088 #nmol−1· h Maximum activity of PGC1a
    Kpg1 = 0.046630145542 #N/A Michaelis-Menten constant for phosphorylation of PGC1a by AMPK
    Kpg2 = 12.3526351747 #N/A Michaelis-Menten constant for deacetylation of PGC1a by SIRT1

    ##for Pulse parameters
    tc1 = 4.38900149 #h Timing of the first AMPK pulse
    tc2 = 15.75 #h Timing of the second AMPK pulse
    tc3 = 18.875 #h Time of maximal nuclear PGC1a abundance
    Td1 = 2.25 #h Duration of the first AMPK pulse
    Td2 = 1.5 #h Duration of the first AMPK pulse
    Td3 = 15.25 #h Duration of the nuclear PGC1a presence
    Tr1 = 2.6 #h Rise time of the first AMPK pulse
    Tr2 = 1.8 #h Rise time of the second AMPK pulse
    Tr3 = 0.5 #h Rise time of nuclear PGC1a
    amp1 = 6.0 #N/A Amplitude of the first AMPK pulse
    amp2 = 0.9778008 #N/A Amplitude of the second AMPK pulse
    amp3 = 0.803062 #nmol−1· l Amplitude of the nuclear PGC1a abundance pulse

    ##paramaters for Chronotherapy timings
    tc4 = 13.664 #Timing of the agonist pulse
    Td4 = 2.83718 #Duration of the agonist pulse
    Tr4 = 1.86794 #Rise time of the agonist pulse
    amp4 = 0.465852 #Amplitude of the agonist pulse

    ##Miscellaneous constants used to describe perturbations
    Csirt = 1 # SIRT1 KO = 0 ; LKB1 KO = 1 ; HFD = 1 ; fasting = 1
    Campk = 1 # SIRT1 KO = 1 ; LKB1 KO = 0.0375 ; HFD = 0.05 ; fasting = 0.05
    offs = 0.02 # SIRT1 KO = 0.02 ; LKB1 KO = 0.02 ; HFD = 0.02 ; fasting = 2.6
    Cpgc1 = 1


    ##for pulse
    #P1=  1/2*((1+np.tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-0)/2.6)) - (1+np.tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-2.25)/2.6)) + (1+np.tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-24)/2.6)) - (1+np.tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-26.25)/2.6)))
    #P2=  1/2*((1+np.tanh((((t-15)- 24*np.absolute((t-15)/24))-0)/1.8)) - (1+np.tanh((((t-15)- 24*np.absolute((t-15)/24))-1.5)/1.8)) + (1+np.tanh((((t-15)- 24*np.absolute((t-15)/24))-24)/1.8)) - (1+np.tanh((((t-15)- 24*np.absolute((t-15)/24))-25.5)/1.8)))
    P1=  1/2*((1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-0)/2.6)) - (1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-2.25)/2.6)) + (1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-24)/2.6)) - (1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-26.25)/2.6)))
    P2=  1/2*((1+np.tanh((((t-15)- 24*((t-15)//24))-0)/1.8)) - (1+np.tanh((((t-15)- 24*((t-15)//24))-1.5)/1.8)) + (1+np.tanh((((t-15)- 24*((t-15)//24))-24)/1.8)) - (1+np.tanh((((t-15)- 24*((t-15)//24))-25.5)/1.8)))
    P3=  1/2*((1+np.tanh((((t-11.25)- 24*((t-11.25)//24))-0)/0.5)) - (1+np.tanh((((t-11.25)- 24*((t-11.25)//24))-15.25)/0.5)) + (1+np.tanh((((t-11.25)- 24*((t-11.25)//24))-24)/0.5)) - (1+np.tanh((((t-11.25)- 24*((t-11.25)//24))-39.25)/0.5)))
    P4=  1/2*((1+np.tanh((((t-12.24541)- 24*((t-12.24541)//24))-0)/1.86794)) - (1+np.tanh((((t-12.24541)- 24*((t-12.24541)//24))-2.83718)/1.86794)) + (1+np.tanh((((t-12.24541)- 24*((t-12.24541)//24))-24)/1.86794)) - (1+np.tanh((((t-12.24541)- 24*((t-12.24541)//24))-26.83718)/1.86794)))



    ##calculation of active ampk, sirt1, pgc1a
    Prot_PGC1a= amp3 * P3
    Act_SIRT = (Csirt*Vsirt*NAD)/(Ksirt+NAD)
    Act_AMPK =Campk* (amp1*P1 + amp2* P2) + (1 - Campk) * offs
    Act_PGC1a= (Cpgc1 * Vpg * Act_AMPK * Act_SIRT * Prot_PGC1a)/(1+((Act_AMPK/Kpg1)*(1+(Act_SIRT/Kpg2))))
    agon_rev= amp4 * P4
    Ki_bmal_rev= (Ki_bmal_rev0/(1+agon_rev))
    Ki_cry_rev= (Ki_cry_rev0/(1+agon_rev))

    PP1.append(Act_AMPK)


    #for per mRNA estimation
    dy1= (Vmax_per*((1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))**hill1)/(1+(CB/(Ka_per_cb*(1+Act_SIRT)))**hill1 * (1+(PC/Ki_per_pc)**hill2))))-dm_per*per

    ##for cry mRNA estimation
    dy2= (Vmax_cry*((1+fold_cry*(CB/(Ka_cry_cb*(1+Act_SIRT)))**hill3)/
                    (1+(CB/(Ka_cry_cb*(1+Act_SIRT)))**hill3 * (1+(PC/Ki_cry_pc)**hill4) * (1+(Prot_rev/Ki_cry_rev)**hill5))))-dm_cry*cry


    ##for rev mRNA estimation
    dy3= (Vmax_rev*((1+fold_rev*(CB/(Ka_rev_cb*(1+Act_SIRT)))**hill6)/(1+(CB/(Ka_rev_cb*(1+Act_SIRT)))**hill6 * (1+(PC/Ki_rev_pc)**hill7))))-dm_rev*rev

    ##for ror mRNA estimation
    dy4= (Vmax_ror*((1+fold_ror*(CB/(Ka_ror_cb*(1+Act_SIRT)))**hill8)/(1+(CB/(Ka_ror_cb*(1+Act_SIRT)))**hill8 * (1+(PC/Ki_ror_pc)**hill9))))-dm_ror*ror

    ##for bmal mRNA estimation
    dy5= (Vmax_bmal*((1+fold_bmal*(1+Act_PGC1a)*(Prot_ror/Ka_bmal_ror)**hill10)/(1+ (Prot_rev/Ki_bmal_rev)**hill11 + (Prot_ror/Ka_bmal_ror)**hill10)))-dm_bmal*bmal

    ##for Per protein estimation
    dy6= Kp_per*per -((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - dp_per*(1+ (m_per_sirt * Act_SIRT) + m_per_ampk * Act_AMPK) * Prot_per

    ##for Cry protein estimation
    dy7= Kp_cry*cry -((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - dp_cry*(1+ m_cry_ampk * Act_AMPK) * Prot_cry

    ##for Rev protein estimation
    dy8= Kp_rev*rev - dp_rev*Prot_rev

    ##for Ror protein estimation
    dy9= Kp_ror*ror - dp_ror*Prot_ror

    ##for bmal protein estimation
    dy10= Kp_bmal*bmal - ((Kass_cb*Prot_bmal)-Kdiss_cb*CB) - dp_bmal*Prot_bmal

    ##for Per-Cry complex protein estimation
    dy11= ((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - d_pc* PC

    ##for Clock-Bmal1 complex protein estimation
    dy12= (Kass_cb*Prot_bmal - Kdiss_cb*CB) - d_cb* CB

    ##for nampt mRNA estimation
    dy13= (Vmax_nampt*((1+fold_nampt*(CB/(Ka_nampt_cb*(1+Act_SIRT)))**hill12)/(1+(CB/(Ka_nampt_cb*(1+Act_SIRT)))**hill12 * (1+(PC/Ki_nampt_pc)**hill13))))-dm_nampt*nampt

    ##for NAMPT Protein estimation
    dy14= (Kp_nampt*nampt)-((dp_nampt* Prot_nampt)/(1+m_nampt_ampk*Act_AMPK))

    ##for NAD level estimation
    dy15= (Vnad*Prot_nampt*((NAD_tot-NAD)/(Knam+NAD_tot-NAD))) - (d_nad*(NAD-NAD_basal)/(Knad+NAD-NAD_basal))

    ##for dbp mRNA estimation
    dy16= (Vmax_dbp*((1+fold_dbp*(CB/(Ka_dbp_cb*(1+Act_SIRT)))**hill14)/(1+(CB/(Ka_dbp_cb*(1+Act_SIRT)))**hill14 * (1+(PC/Ki_dbp_pc)**hill15))))-dm_dbp*dbp


    return (dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16)


""""[per_mrna,cry_mrna,rev_mrna,ror_mrna,bmal_mrna,Prot_per,  Prot_cry ,      Prot_rev   ,   Prot_ror      ,Prot_bmal,      PC,      CB     ,nampt_mrna,  NAMPT,   NAD+, dbp_mrna]"""
#x = [0.45,   1.35,     0.25,    0.994,   2.35,        0.256,       1.37,     0.0455679554,       1.20633514,     3.27,   6.49    ,    0.3410,      0.45,   0.75,    1.2,  0.154]
x=[1.8260362923896098, 1.1793403022645064, 0.2845456944513496, 0.4866078167898578, 0.06921850482736103, 1.8382769791479783, 0.1958873081880534, 0.16463272427500814, 1.0064165976483035, 0.9488651725940949, 4.632529156858846, 0.1837328621162798, 1.2310606314850925, 1.8161399412996202, 1.6106878781932097, 0.8963281546317611]

ff=900
tspan=[0, ff];
kk= ff*10
t_eval = np.linspace(0, ff, kk)

#t_conver= time_conver(t_eval,98.0)
yd = solve_ivp(fun= odes2, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=x)

'''

#to create diff subplots for per, cry, rev, ror mrna

f, (plt1, plt2, plt3, plt4) = plt.subplots(4) #creating diff subplots

line1 = plt1.plot(t_eval, yd['y'][0], label="per mrna", color="blue") #dy1
#plot time points, at each timepoint
#saying we want every row in the first column to plot each timepoint
line2 = plt2.plot(t_eval, yd['y'][1], label="cry mrna", color="green") #dy2
line3 = plt3.plot(t_eval, yd['y'][2], label="rev mrna", color="red")#dy3
line4 = plt4.plot(t_eval, yd['y'][3], label="ror mrna", color="orange") #dy4



plt4.set_xlabel('Time')
plt1.set_ylabel('per')
plt2.set_ylabel('cry')
plt3.set_ylabel('rev')
plt4.set_ylabel('ror')

'''

'''
#Core Clock Loop
per cry protein complex and clock bmal protein complex
#per/cry complex inhibits clock:bmal1 in the core feedback loop


f, (plt1, plt2) = plt.subplots(2) #creating diff subplots

line1 = plt1.plot(t_eval, yd['y'][10], color="blue") #dy11

line2 = plt2.plot(t_eval, yd['y'][11], color="green") #dy12

plt2.set_xlabel('Time')
plt1.set_ylabel('PER CRY complex')
plt2.set_ylabel('CLOCK BMAL1 complex')



#per, cry, bmal1 mrna in one plot

plt.plot(t_eval,yd['y'][0],label ='per mrna', color = 'green')
plt.plot(t_eval,yd['y'][1], label = 'cry mrna', color = 'blue')
plt.plot(t_eval,yd['y'][4], label = 'bmal1 mrna', color = 'orange')
plt.ylabel('Concentration')
plt.xlabel('Time')
plt.legend()


'''
'''

#rev-erb, ror feedback on bmal1 (antiphase rhythm)


plt.plot(t_eval, yd['y'][7], 'r', label="REV-ERB protein")

plt.plot(t_eval, yd['y'][8], 'g', label="ROR protein")

plt.plot(t_eval, yd['y'][4], 'b', label="BMAL1 mRNA")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.title("REV-ERB and ROR feedback on BMAL1 transcription")
plt.legend()
plt.show()


'''
'''

#Rev-erb inhibition on Cry1 and Bmal1 protein


f, (plt1, plt2, plt3) = plt.subplots(3) #creating diff subplots

line1 = plt1.plot(t_eval, yd['y'][6], label="Cry1", color="blue") #dy7

line2 = plt2.plot(t_eval, yd['y'][9], label="Bmal1", color="green") #dy10

line3 = plt3.plot(t_eval, yd['y'][7], label="REV-ERB", color="red") #dy8

plt3.set_xlabel('Time')
plt1.set_ylabel('Cry1')
plt2.set_ylabel('Bmal1')
plt3.set_ylabel('REV ERB')

# how to plot AMPK and SIRT2 relationship


hh = []   # AMPK activity
sirt = [] # SIRT1 activity

for i, t in enumerate(t_eval):

    #First AMPK pulse
    P1 = 0.5*((1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-0)/2.6))
              - (1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-2.25)/2.6))
              + (1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-24)/2.6))
              - (1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-26.25)/2.6)))

    #Second AMPK pulse
    P2_val = 0.5*((1+np.tanh(((t-15)-24*((t-15)//24)-0)/1.8))
                  - (1+np.tanh(((t-15)-24*((t-15)//24)-1.5)/1.8))
                  + (1+np.tanh(((t-15)-24*((t-15)//24)-24)/1.8))
                  - (1+np.tanh(((t-15)-24*((t-15)//24)-25.5)/1.8)))

    #Active AMPK signal
    Campk = 1
    amp1, amp2, offs = 6.0, 0.9778008, 0.02
    Act_AMPK = Campk*(amp1*P1 + amp2*P2_val) + (1-Campk)*offs
    hh.append(Act_AMPK)

    #Active SIRT1 signal(depends on NAD)
    NAD = yd['y'][14][i]  # NAD
    Csirt, Vsirt, Ksirt = 1.0, 1.0, 1.0  #parameter values
    Act_SIRT = (Csirt*Vsirt*NAD)/(Ksirt+NAD)
    sirt.append(Act_SIRT)

#line plot with AMPK and SIRT activity
plt.figure(figsize=(8,5))
plt.plot(t_eval, hh, 'r', label="AMPK activity")
plt.plot(t_eval, sirt, 'g', label="SIRT1 activity")
plt.xlabel("Time (h)", fontweight='bold')
plt.ylabel("Activity", fontweight='bold')
plt.legend()
plt.show()


'''
'''
#metabolic module - AMPK, SIRT1, NAD, NAMPT


nampt_mrna = yd['y'][12]
nampt_prot = yd['y'][13]
nad = yd['y'][14]

# Compute SIRT1 activity
Csirt, Vsirt, Ksirt = 1.0, 1.0, 1.0   #parameter values
sirt = (Csirt*Vsirt*nad)/(Ksirt+nad)

#AMPK activity using pulses
hh = []
for t in t_eval:
    P1 = 0.5*((1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-0)/2.6))
              - (1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-2.25)/2.6))
              + (1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-24)/2.6))
              - (1+np.tanh(((t-3.26400149)-24*((t-3.26400149)//24)-26.25)/2.6)))
    P2_val = 0.5*((1+np.tanh(((t-15)-24*((t-15)//24)-0)/1.8))
                  - (1+np.tanh(((t-15)-24*((t-15)//24)-1.5)/1.8))
                  + (1+np.tanh(((t-15)-24*((t-15)//24)-24)/1.8))
                  - (1+np.tanh(((t-15)-24*((t-15)//24)-25.5)/1.8)))
    amp1, amp2, offs, Campk = 6.0, 0.9778008, 0.02, 1
    Act_AMPK = Campk*(amp1*P1 + amp2*P2_val) + (1-Campk)*offs
    hh.append(Act_AMPK)


plt.plot(t_eval, hh, 'r', label="AMPK activity")
plt.plot(t_eval, sirt, 'g', label="SIRT1 activity")
plt.plot(t_eval, nampt_prot, 'b', label="NAMPT protein")
plt.plot(t_eval, nad, 'm', label="NAD+")
plt.xlabel("Time")
plt.ylabel("Activity")
plt.title("Metabolic Module Dynamics")
plt.legend()
plt.show()

#PGC1 activated by both AMPK and SIRT1

#end of my code


'''
'''
#plt.plot(t_eval,yd['y'][6], 'b')
plt.plot(t_eval,yd['y'][0], 'g')
plt.xlabel("Time", fontweight='bold')
plt.ylabel("Conc",fontweight='bold')
plt.legend()

#plt.xticks(tconver2, fontweight='bold')
plt.show()




plt.plot(t_eval,yd['y'][13], 'g')
plt.plot(t_eval,yd['y'][14], 'b')
plt.xlabel("Time", fontweight='bold')
plt.ylabel("Conc",fontweight='bold')
plt.legend()



#


#print (type(abs(-5.5)))
hh=[]
P2=[]
##for active ampk
for t in t_eval:
    tc1 = 4.38900149 #h Timing of the first AMPK pulse
    tc2 = 15.75 #h Timing of the second AMPK pulse
    tc3 = 18.875 #h Time of maximal nuclear PGC1a abundance
    Td1 = 2.25 #h Duration of the first AMPK pulse
    Td2 = 1.5 #h Duration of the first AMPK pulse
    Td3 = 15.25 #h Duration of the nuclear PGC1a presence
    Tr1 = 2.6 #h Rise time of the first AMPK pulse
    Tr2 = 1.8 #h Rise time of the second AMPK pulse
    Tr3 = 0.5 #h Rise time of nuclear PGC1a
    amp1 = 6.0 #N/A Amplitude of the first AMPK pulse
    amp2 = 0.9778008 #N/A Amplitude of the second AMPK pulse
    amp3 = 0.803062 #nmol−1· l Amplitude of the nuclear PGC1a abundance pulse

    Campk = 1
    amp1 = 6.0
    amp2 = 0.9778008
    offs = 0.02

    P1= 1/2*((1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-0)/2.6)) - (1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-2.25)/2.6)) + (1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-24)/2.6)) - (1+np.tanh((((t-3.26400149)- 24*((t-3.26400149)//24))-26.25)/2.6)))
    P2.append(1/2*((1+np.tanh((((t-15)- 24*((t-15)//24))-0)/1.8)) - (1+np.tanh((((t-15)- 24*((t-15)//24))-1.5)/1.8)) + (1+np.tanh((((t-15)- 24*((t-15)//24))-24)/1.8)) - (1+np.tanh((((t-15)- 24*((t-15)//24))-25.5)/1.8))))

    #hh.append(Campk* (amp1*P1 + amp2* P2) + (1 - Campk) * offs)
    hh.append(P1)

plt.plot(t_eval[5016:5417],hh[5016:5417], 'r', label= "Active AMPK",markersize=3.5)
plt.plot(t_eval[5016:5417],yd['y'][13][5016:5417], 'g', label= "NAMPT",markersize=3.5)
plt.plot(t_eval[5016:5417],yd['y'][14][5016:5417], 'b',label= "NAD+",markersize=3.5)
plt.axvspan(501.656,511.157, color='gray')
plt.axvspan(523.858,535.159, color='gray')
plt.xlabel("Time", fontweight='bold')
plt.ylabel("Conc",fontweight='bold')
plt.legend()
#plt.xticks(tconver2, fontweight='bold')
plt.show()




fig, ax = plt.subplots(1,1, figsize=[11, 9])




#plt.subplot(3,2,1)
ax[0,0].plot(t_eval[5016:5417],yd['y'][13][5016:5417], 'g', label= "NAMPT")
ax[0,0].axvspan(5016,5111, color='gray')
ax[0,0].axvspan(5238,5351, color='gray')
ax[0,0].legend()







'''

