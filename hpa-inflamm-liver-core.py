import math as mt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
%matplotlib
import matplotlib.pyplot as plt
import random
import pandas as pd

def liver_core_Hpa_infla(t,y):

# the HPA axis and inflamaton

    CRH=y[0];
    ACTH=y[1];
    StARp=y[2];
    CORT=y[3];
    Dex1 =y[4];
    Dex2 =y[5];

    CORTp =y[6];

    GR_mrna =y[7];
    GR_prot =y[8]; #total protein
    GR_cyt =y[9];
    GR =y[10]; #nucleus

    LPS=y[11];
    Phg=y[12];
    Phg1= y[13];
    TGF=y[14];
    TNF=y[15];
    IL10=y[16];
    IL6=y[17];

    # From here liver core model starts

    per= y[18]
    cry= y[19]
    rev= y[20]
    ror= y[21]
    bmal= y[22]
    Prot_per= y[23]
    Prot_cry= y[24]
    Prot_rev= y[25]
    Prot_ror= y[26]
    Prot_bmal= y[27]
    PC= y[28]
    CB= y[29]
    nampt= y[30]
    Prot_nampt= y[31]
    NAD= y[32]
    dbp= y[33]


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
    dx=1; #baseline dose value
    dex=0;      #Comment this to run IC50 test


    d=0; #shift/delay parameter (used to offset dosing times)

    if inputF == "no":
        if t >= 98.25+d and t < 98.5+d or t >= 122.0+d and t < 122.25+d or t >=146.0+d and t < 146.25+d or t >= 170.0+d and t < 170.25+d or t >= 194.0+d and t < 194.25+d or t >= 218.0+d and t < 218.25+d or t >= 242.0+d and t < 242.25+d or t >= 266.0+d and t < 266.25+d:
            dex=dx*10;      #Comment this to run IC50 test
        #if t >= 98.0+d and t < 98.25+d or t >= 122.0+d and t < 122.25+d or t >=146.0+d and t < 146.25+d or t >= 170.0+d and t < 170.25+d or t >= 194.0+d and t < 194.25+d or t >= 218.0+d and t < 218.25+d or t >= 242.0+d and t < 242.25+d or t >= 266.0+d and t < 266.25+d: dex=dx*10;      #Comment this to run IC50 test

    #here if d is 0, than the pulses will happen at 1 am everyday at that interval
    #if d is 1, then it will hapen at 2 am, etc.

    Dex=dex           #*((t>696.00+d) & (t<696.25+d));



    '''%%%%%%%%%%%%%%%% HPA-Axis %%%%%%%%%%%%%'''


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



    # parameter values from @ Mavroudis et.al. 2012
    kc = 0.009  #1/hr    coupling parameter,  coupling strength between central and peripheral circadian systems
    kac = 5     #acetylation rate
    #IC50rm = 26.2 # synthesis drop to half

    # Circadian drive from SCN clock
    omega=2*mt.pi/24;
    cir=2*(1+ mt.cos(omega*(t)));

    # Circadian cdrive from adrenal peripheral clock
    cir2=(4-cir)*(Ki_GR/(Ki_GR+GR));


    # Corticotrophic releasing hormone (CRH)
    dy1 = (K_strs*(1+st)*(Ki_GR1**n1/(Ki_GR1**n1+GR**n1))*(cir)*(1+V_crh*TNF)-k_crh*CRH)

    # Adrenocorticotropic hormone (ACTH)
    dy2 = (V_acth1*CRH*(Ki_GR1**n1/(Ki_GR1**n1+GR**n1))*(1+V_acth2*((TNF)**n2/(Km_tnf3**n2+ TNF**n2)))-k_acth*ACTH)

    # StAR proetin
    dy3 =s0*(V_starp1*(ACTH*cir2)*(1+V_starp2*(TNF**n2/(Km_tnf3**n2+TNF**n2))*(Kil6/(Kil6+IL6)))  -k_starp*StARp)

    # Cortisol
    dy4 = (V_cort*StARp-k_cort*CORT)

    # Dexamethasone kinetics
    dy5 = Dex-k_dex1*Dex1
    dy6 = V_dex2*Dex1-k_dex2*Dex2

    # Peripheral cortisol
    dy7 =(1/(tc))*(CORT+(fac*Dex2)-CORTp)

    # GR mRNA
    dy8 = V_GRmrna*(1 - ((GR/(kac*CB))**n3/((km_GR1)**n3+(GR/(kac*CB))**n3)))- k_GRmrna*GR_mrna; #added CB modification
    #GR binding to GRE is inhibited by the action of the CLOCK/BMAL1 heterodimer that acetylizes GR nuclear complex in its hinge region and blocks its forward action
    #dy(8) = V_GRmrna*(1 - (GR**n3/((km_GR1)**n3+GR**n3)))- k_GRmrna*GR_mrna;

    # GR protein
    dy9 = V_GRprot*GR_mrna + fr*kre*GR- kon*(CORTp)*GR_prot - k_GRprot*GR_prot

    # GR-cortisol complex in cytosol
    dy10 = kon*(CORTp)*GR_prot- krt*GR_cyt*(ki_tnf**n4/(TNF**n4+ki_tnf**n4))

    # Nuclear GR
    dy11 = krt*GR_cyt*(ki_tnf**n4/(TNF**n4+ki_tnf**n4)) -kre*GR



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
    cir3=2*(1+ mt.cos(omega*(t)+0.5)) #rhythym in peripheral phagocytic cells

    # LPS
    dy12 = 1e-7*(1+It)-k_lps*LPS*Phg;

    # Phagocytes

    dy13 = (V_phg1*(((1*cir3)+ (V_phg2* TNF/(Km_tnf1+TNF)))* (Ki_tgf1/(Ki_tgf1+ TGF))* (Ki_IL10/(Ki_IL10+IL10)))*LPS- k_phg*Phg);

    dy14 = (1.0/(tp))*(Phg-Phg1);

    # Transformig grwth factor (TGF)
    dy15 =  (V_tgf1* Phg + (V_tgf2*GR**n5/(Km_GR2**n5+GR**n5))- k_tgf*TGF);

    # Tumour necrosis factor (TNF)
    dy16 = tnf_b+(Phg/(Km_phg1+Phg))*(V_tnf1+V_tnf2*(TNF/(Km_tnf2+TNF)))*(Ki_tgf2**n6/(Ki_tgf2**n6+TGF**n6))*(1-(IL6**n7/(Km_il61**n7+IL6**n7)))-k_tnf*(TNF);

    # Interleukin 10 (IL10)
    dy17 = ((V_il101*Phg1**n8/(Phg1**n8+Km_phg2**n8))+ (V_il102* TGF**n9/(Km_tgf3**n9+ TGF**n9)) - k_il10*IL10);

    # Interleukin 6 (IL6)
    dy18 = il6_b+(Phg1/(Km_phg3+Phg1))*(V_il61+V_il62*(TNF+IL6)**n10/(Km_il6tnf**n10+(TNF+IL6)**n10))*(Ki_il102/(Ki_il102+(IL10)))-k_il6*IL6;



    #parameters for mRNA and protein degradatation rates
    dm_bmal = 0.827333442085  #Bmal mRNA degradation rate
    dm_cry = 0.319848706181   #Cry mRNA degradation rate
    dm_dbp = 0.384486532062   #Dbp mRNA degradation rate
    dm_nampt = 0.814311309051 #Nampt mRNA degradation rate
    dm_per = 0.323114598647   #Per mRNA degradation rate
    dm_rev = 4.0479072846     #Rev-Erb mRNA degradation rate
    dm_ror = 0.246575760727   #Ror mRNA degradation rate
    dp_bmal = 0.186413466797  #BMAL protein degradation rate
    dp_cry = 0.599026119971   #CRY protein degradation rate
    dp_nampt = 49.8841023982  #NAMPT protein degradation rate
    dp_per = 10.9446515392    #PER protein degradation rate
    dp_rev = 0.281564063057   #REV-ERB protein degradation rate
    dp_ror = 0.0340112196281  #ROR degradation rate
    d_cb =0.197714012552      #CLOCK-BMAL complex degradation rate
    d_pc = 0.609290138329     #PER-CRY complex degradation rate

    #parameters for Complexation kinetic rates
    Kass_cb = 0.0162923358655    #nmol−1· l · h−1 CLOCK-BMAL association rate
    Kass_pc = 12.302005485       #nmol−1· l · h−1 PER-CRY association rate
    Kdiss_cb = 0.00613502224231  #h−1 CLOCK-BMAL dissociation rate
    Kdiss_pc = 0.0365867175408   #h−1 PER-CRY dissociation rate


    #parmeters for Maximal transcription rates
    Vmax_bmal = 0.0916862770193 #Bmal maximal transcription rate
    Vmax_cry = 0.702216664807 #Cry maximal transcription rate
    Vmax_dbp = 0.0802009609453 #Dbp maximal transcription rate
    Vmax_nampt = 3.49035201578 #Nampt maximal transcription rate
    Vmax_per = 0.730201742662 #Per maximal transcription rate
    Vmax_rev = 1.12297601784 #Rev-Erb maximal transcription rate
    Vmax_ror = 6.9843472736 #Ror maximal transcription rate

    #parameters for Activation ratios
    fold_bmal = 15.9258093373 #activation ratio of Bmal by ROR
    fold_cry = 1.1604489571 #activation ratio of Cry by CLOCK-BMAL
    fold_dbp = 400.0 #activation ratio of Dbp by CLOCK-BMAL
    fold_nampt = 1.57880681573 #activation ratio of Nampt by CLOCK-BMAL
    fold_per = 12.977351123 #activation ratio of Per by CLOCK-BMAL
    fold_rev = 73.2342431701 #activation ratio of Rev-Erb by CLOCK-BMAL
    fold_ror = 335.923333883 #activation ratio of Ror by CLOCK-BMAL

    #parameters for
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

    #parameters for Hill coefficients
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

    #parameters for Translation rates
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

    #parameters for NAD kinetics, Sirt1 and PGC1a activity
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


    #for Pulse parameters
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

    #paramaters for Chronotherapy timings
    tc4 = 13.664 #Timing of the agonist pulse
    Td4 = 2.83718 #Duration of the agonist pulse
    Tr4 = 1.86794 #Rise time of the agonist pulse
    amp4 = 0.465852 #Amplitude of the agonist pulse

    #Miscellaneous constants used to describe perturbations
    Csirt = 1 # SIRT1 KO = 0 ; LKB1 KO = 1 ; HFD = 1 ; fasting = 1
    Campk = 1 # SIRT1 KO = 1 ; LKB1 KO = 0.0375 ; HFD = 0.05 ; fasting = 0.05
    offs = 0.02 # SIRT1 KO = 0.02 ; LKB1 KO = 0.02 ; HFD = 0.02 ; fasting = 2.6
    Cpgc1 = 1


    #for pulse
    #P1=  1/2*((1+tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-0)/2.6)) - (1+tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-2.25)/2.6)) + (1+tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-24)/2.6)) - (1+tanh((((t-3.26400149)- 24*np.absolute((t-3.26400149)/24))-26.25)/2.6)))
    #P2=  1/2*((1+tanh((((t-15)- 24*np.absolute((t-15)/24))-0)/1.8)) - (1+tanh((((t-15)- 24*np.absolute((t-15)/24))-1.5)/1.8)) + (1+tanh((((t-15)- 24*np.absolute((t-15)/24))-24)/1.8)) - (1+tanh((((t-15)- 24*np.absolute((t-15)/24))-25.5)/1.8)))
    P1 =  1 #1/2*((1+tanh((((t-3.26400149)- 24*floor((t-3.26400149)/24))-0)/2.6)) - (1+tanh((((t-3.26400149)- 24*floor((t-3.26400149)/24))-2.25)/2.6)) + (1+tanh((((t-3.26400149)- 24*floor((t-3.26400149)/24))-24)/2.6)) - (1+tanh((((t-3.26400149)- 24*floor((t-3.26400149)/24))-26.25)/2.6))) ;
    P2 =  1 #1/2*((1+tanh((((t-15)- 24*floor((t-15)/24))-0)/1.8)) - (1+tanh((((t-15)- 24*floor((t-15)/24))-1.5)/1.8)) + (1+tanh((((t-15)- 24*floor((t-15)/24))-24)/1.8)) - (1+tanh((((t-15)- 24*floor((t-15)/24))-25.5)/1.8))) ;
    P3 =  1 #1/2*((1+tanh((((t-11.25)- 24*floor((t-11.25)/24))-0)/0.5)) - (1+tanh((((t-11.25)- 24*floor((t-11.25)/24))-15.25)/0.5)) + (1+tanh((((t-11.25)- 24*floor((t-11.25)/24))-24)/0.5)) - (1+tanh((((t-11.25)- 24*floor((t-11.25)/24))-39.25)/0.5))) ;
    P4 =  1  #1/2*((1+tanh((((t-12.24541)- 24*floor((t-12.24541)/24))-0)/1.86794)) - (1+tanh((((t-12.24541)- 24*floor((t-12.24541)/24))-2.83718)/1.86794)) + (1+tanh((((t-12.24541)- 24*floor((t-12.24541)/24))-24)/1.86794)) - (1+tanh((((t-12.24541)- 24*floor((t-12.24541)/24))-26.83718)/1.86794))) ;
    #changing all pulses to 1, will be static

    fFA = 1 #free fatty acids
    fAA = 1 #free amino acids
    Glup = 1 #glucose uptake

    # calculation of active ampk, sirt1, pgc1a
    Prot_PGC1a= amp3 * P3 #directly scaled
    Act_SIRT = (Csirt*Vsirt*NAD)/(Ksirt+NAD)  # MM type equation

    #Act_AMPK =Campk* (amp1*P1 + amp2*P2) + (1 - Campk) * offs ;

    Act_AMPK = (-0.056739 * fFA + 0.018007 * (fFA ** 2) + 0.00053297 * (fFA ** 3) + 0.78872 * fAA
        - 0.03337 * fAA * fFA - 0.010077 * fAA * (fFA ** 2) + 0.0010072 * fAA * (fFA ** 3) - 0.65065 * (fAA ** 2)
        + 0.033309 * (fAA ** 2) * fFA + 0.00098831 * (fAA ** 2) * (fFA ** 2) + 0.21443 * (fAA ** 3) - 0.0059495 * (fAA ** 3) * fFA
        - 0.1838 * Glup + 0.12119 * Glup * fFA - 0.016247 * Glup * (fFA ** 2)  + 0.0011439 * Glup * (fFA ** 3) + 0.31752 * Glup * fAA
        - 0.038471 * Glup * fAA * fFA + 2.1693e-05 * Glup * fAA * (fFA ** 2) - 0.19239 * Glup * (fAA ** 2) + 0.0069602 * Glup * (fAA ** 2) * fFA
        + 0.02598 * Glup * (fAA ** 3) + 0.12704 * (Glup ** 2) - 0.019491 * (Glup ** 2) * fFA + 0.0015665 * (Glup ** 2) * (fFA ** 2) + 0.029975 * (Glup ** 2) * fAA
        - 5.8394e-06 * (Glup ** 2) * fAA * fFA - 0.00033581 * (Glup ** 2) * (fAA ** 2) - 0.058545 * (Glup ** 3) + 0.0019572 * (Glup ** 3) * fFA
        - 0.0032067 * (Glup ** 3) * fAA - 0.086452 * 1 + 0.0072993 * (Glup ** 4) - 0.023548 * (fAA ** 4) - 0.00044542 * (fFA ** 4))


    Act_PGC1a= (Cpgc1 * Vpg * Act_AMPK * Act_SIRT * Prot_PGC1a)/(1+((Act_AMPK/Kpg1)*(1+(Act_SIRT/Kpg2))))
    agon_rev= amp4 * P4
    Ki_bmal_rev= (Ki_bmal_rev0/(1+agon_rev))
    Ki_cry_rev= (Ki_cry_rev0/(1+agon_rev))

    cir4= (1+ mt.cos(omega*(t)+7.75))        # 1;

    #for per mRNA estimation
    dy19 = 1.5* (Vmax_per*((1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))**hill1)/(1+(CB/(Ka_per_cb*(1+Act_SIRT)))**hill1 * (1+(PC/Ki_per_pc)**hill2))))-dm_per*per + 1.5*kc*(GR/(kac*CB*100))
    #dy(19)= (Vmax_per*((1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))**hill1)/(1+(CB/(Ka_per_cb*(1+Act_SIRT)))**hill1 * (1+(PC/Ki_per_pc)**hill2))))-dm_per*per;

    #for cry mRNA estimation
    dy20 = (Vmax_cry*((1+fold_cry*(CB/(Ka_cry_cb*(1+Act_SIRT)))**hill3)/(1+(CB/(Ka_cry_cb*(1+Act_SIRT)))**hill3 * (1+(PC/Ki_cry_pc)**hill4) * (1+(Prot_rev/Ki_cry_rev)**hill5))))-dm_cry*cry

    #for rev mRNA estimation
    dy21 = (Vmax_rev*((1+fold_rev*(CB/(Ka_rev_cb*(1+Act_SIRT)))**hill6)/(1+(CB/(Ka_rev_cb*(1+Act_SIRT)))**hill6 * (1+(PC/Ki_rev_pc)**hill7))))-dm_rev*rev

    #For ror mRNA estimation
    dy22 = (Vmax_ror*((1+fold_ror*(CB/(Ka_ror_cb*(1+Act_SIRT)))**hill8)/(1+(CB/(Ka_ror_cb*(1+Act_SIRT)))**hill8 * (1+(PC/Ki_ror_pc)**hill9))))-dm_ror*ror

    #for bmal mRNA estimation
    dy23 = (0.45*Vmax_bmal*((1+fold_bmal*(1+Act_PGC1a)*(Prot_ror/Ka_bmal_ror)**hill10)/(1+ (Prot_rev/Ki_bmal_rev)**hill11 + (Prot_ror/Ka_bmal_ror)**hill10)))*cir4-dm_bmal*bmal

    #for Per protein estimation
    dy24 = Kp_per*per -((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - dp_per*(1+ (m_per_sirt * Act_SIRT) + m_per_ampk * Act_AMPK) * Prot_per

    #for Cry protein estimation
    dy25 = Kp_cry*cry -((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - dp_cry*(1+ m_cry_ampk * Act_AMPK) * Prot_cry

    #for Rev protein estimation
    dy26 = Kp_rev*rev - dp_rev*Prot_rev

    #for Rev protein estimation
    dy27 = Kp_ror*ror - dp_ror*Prot_ror

    #for bmal protein estimation
    dy28 = Kp_bmal*bmal - ((Kass_cb*Prot_bmal)-Kdiss_cb*CB) - dp_bmal*Prot_bmal

    #for Per-Cry complex protein estimation
    dy29 = ((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - d_pc* PC

    #for Clock-Bmal1 complex protein estimation
    dy30 = (Kass_cb*Prot_bmal - Kdiss_cb*CB) - d_cb* CB

    #for nampt mRNA estimation
    dy31 = (Vmax_nampt*((1+fold_nampt*(CB/(Ka_nampt_cb*(1+Act_SIRT)))**hill12)/(1+(CB/(Ka_nampt_cb*(1+Act_SIRT)))**hill12 * (1+(PC/Ki_nampt_pc)**hill13))))-dm_nampt*nampt

    #for NAMPT Protein estimation
    dy32 = (Kp_nampt*nampt)-((dp_nampt* Prot_nampt)/(1+m_nampt_ampk*Act_AMPK))

    #for NAD level estimation
    dy33 = (Vnad*Prot_nampt*((NAD_tot-NAD)/(Knam+NAD_tot-NAD))) - (d_nad*(NAD-NAD_basal)/(Knad+NAD-NAD_basal))

    #for dbp mRNA estimation
    dy34 = (Vmax_dbp*((1+fold_dbp*(CB/(Ka_dbp_cb*(1+Act_SIRT)))**hill14)/(1+(CB/(Ka_dbp_cb*(1+Act_SIRT)))**hill14 * (1+(PC/Ki_dbp_pc)**hill15))))-dm_dbp*dbp

    return (dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16,dy17,dy18,dy19,dy20,dy21,dy22,dy23,dy24,dy25,dy26,dy27,dy28,dy29,dy30,dy31,dy32,dy33,dy34)



p=0.0; #LPS injection
r=1.0; # Change this to vary GR sensitivity to cortisol
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


#inputF = 'yes' #no dexamethasone for ic50 test
inputF = 'no' #dexamethasone for ic50 test


y0=[24.95694763498119,4.9123730894505435,0.9938828562096904,7.393415730911929,0.0,0.0,7.636200249054682,18.001529048278428,350.07408500237455,19.87971898349494,28.182728997762037,6.65145308687877e-05,41.97935885594562,50.96107608767422,0.2872674163673849,1.0053860612262895,3.5041839210619,1.072968514441532,1.8260362923896098, 1.1793403022645064, 0.2845456944513496, 0.4866078167898578, 0.06921850482736103, 1.8382769791479783, 0.1958873081880534, 0.16463272427500814, 1.0064165976483035, 0.9488651725940949, 4.632529156858846, 0.1837328621162798, 1.2310606314850925, 1.8161399412996202, 1.6106878781932097, 0.8963281546317611]

yd = solve_ivp(fun= liver_core_Hpa_infla, method = 'BDF',t_span=tspan, t_eval= t_eval , y0=y0 , rtol = 1e-10,atol = 1e-10)

''' %%%

#metabolic module - AMPK, SIRT1, NAD, NAMPT
#PGC1 activated by both AMPK and SIRT1

t = yd.t
nampt = yd.y[30, :]
Prot_nampt = yd.y[31, :]
NAD = yd.y[32, :]

#computing activity of SIRT - constants copied from the ODE function
Vsirt = 0.915854846427
Ksirt = 0.75

Act_SIRT = (Vsirt * NAD) / (Ksirt + NAD)

#Act_AMPK is a constant because fFA = fAA = Glup = 1

fFA = 1.0
fAA = 1.0
Glup = 1.0

Act_AMPK = (-0.056739 * fFA + 0.018007 * (fFA ** 2) + 0.00053297 * (fFA ** 3) + 0.78872 * fAA
    - 0.03337 * fAA * fFA - 0.010077 * fAA * (fFA ** 2) + 0.0010072 * fAA * (fFA ** 3) - 0.65065 * (fAA ** 2)
    + 0.033309 * (fAA ** 2) * fFA + 0.00098831 * (fAA ** 2) * (fFA ** 2) + 0.21443 * (fAA ** 3) - 0.0059495 * (fAA ** 3) * fFA
    - 0.1838 * Glup + 0.12119 * Glup * fFA - 0.016247 * Glup * (fFA ** 2)  + 0.0011439 * Glup * (fFA ** 3) + 0.31752 * Glup * fAA
    - 0.038471 * Glup * fAA * fFA + 2.1693e-05 * Glup * fAA * (fFA ** 2) - 0.19239 * Glup * (fAA ** 2) + 0.0069602 * Glup * (fAA ** 2) * fFA
    + 0.02598 * Glup * (fAA ** 3) + 0.12704 * (Glup ** 2) - 0.019491 * (Glup ** 2) * fFA + 0.0015665 * (Glup ** 2) * (fFA ** 2) + 0.029975 * (Glup ** 2) * fAA
    - 5.8394e-06 * (Glup ** 2) * fAA * fFA - 0.00033581 * (Glup ** 2) * (fAA ** 2) - 0.058545 * (Glup ** 3) + 0.0019572 * (Glup ** 3) * fFA
    - 0.0032067 * (Glup ** 3) * fAA - 0.086452 * 1 + 0.0072993 * (Glup ** 4) - 0.023548 * (fAA ** 4) - 0.00044542 * (fFA ** 4))

#PGC1a activity
Vpg = 24.06372088
Kpg1 = 0.046630145542
Kpg2 = 12.3526351747
Cpgc1 = 1.0
amp3 = 0.803062
P3 = 1.0  # in your current code
Prot_PGC1a = amp3 * P3

Act_PGC1a = (Cpgc1 * Vpg * Act_AMPK * Act_SIRT * Prot_PGC1a) / (1 + ((Act_AMPK / Kpg1) * (1 + (Act_SIRT / Kpg2))))


plt.plot(t, NAD, label='NAD')
plt.plot(t, Prot_nampt, label='NAMPT protein')
plt.plot(t, Act_SIRT, label='SIRT1 activity')
plt.hlines(Act_AMPK, t[0], t[-1], colors='r', linestyles='--', label='AMPK activity (const)')
plt.hlines(Act_PGC1a, t[0], t[-1], colors='b', linestyles='--', label='PGC1a activity (const)')
plt.xlabel('Time (h)')
plt.ylabel('Value')
plt.legend()
plt.show()



t = yd.t                          # all 9000 time points in an array
timeframe = t[4000:]              # time points from 4000 onwards
y_frame = yd.y[:, 4000:]          # shape (23, len(timeframe))

#Peaks
peaks_indices = [find_peaks(y_frame[i])[0].tolist()   # shape[0] refers to the number of rows, [0] to only take index
                 for i in range(y_frame.shape[0])]
peaks_values  = [y_frame[i, find_peaks(y_frame[i])[0]].tolist()
                 for i in range(y_frame.shape[0])]

#Troughs
troughs_indices = [find_peaks(-y_frame[i])[0].tolist()  # we use negative
                   for i in range(y_frame.shape[0])]
troughs_values  = [y_frame[i, find_peaks(-y_frame[i])[0]].tolist()
                   for i in range(y_frame.shape[0])]


peaks_times = [timeframe[idx].tolist() for idx in peaks_indices]
troughs_times = [timeframe[idx].tolist() for idx in troughs_indices] # true time from 0 to 900

print(yd.y[6, 4171])


# to get a table for a certain value, just change the number

Peak_Values = peaks_values[6]
Peak_Times = peaks_times[6]
peak_df = pd.DataFrame({ "Peak Value": Peak_Values, "Peak Time": Peak_Times })
#print(peak_df)

Trough_Values = troughs_values[6]
Trough_Times = troughs_times[6]
trough_df = pd.DataFrame({ "Trough Value": Trough_Values, "Trough Time": Trough_Times })
#print(trough_df)

'''
