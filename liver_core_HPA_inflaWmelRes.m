function df= liver_core_HPA_inflaWmelRes(t,y)
%% the HPA axis and inflamaton


%global n  Ki dex nj R Ki1  p r s0  u Q inputF d bcz
global  n  Ki dex nj R Ki1  p r s0  u Q p1 inputF d bcz fFA fAA Glup
 
CRH=y(1);
ACTH=y(2);
StARp=y(3);
CORT=y(4);
Dex1 =y(5);
Dex2 =y(6);

CORTp =y(7);
GR_mrna =y(8);
GR_prot =y(9);
GR_cyt =y(10);
GR =y(11);

LPS=y(12);
Phg=y(13);
Phg1= y(14);
TGF=y(15);
TNF=y(16);
IL10=y(17);
IL6=y(18);

%%From here meal reasponse varaiables starts
Qst1=y(19);
Qst2=y(20);
Glu_gut=y(21);
Pro_gut=y(22);
Pro_Int=y(23);
Fat_gut=y(24);
Fat_Int=y(25);
InsV=y(26);
InsL=y(27);
InsP=y(28);
Cb_Ala_L=y(29);
Cb_FFA_L=y(30);
Ca_Glup	=y(31);
Ca_Glu_OT=y(32);
Ca_Alap=y(33);
Ca_FFAp=y(34);

%%From here liver core model strats

%print (t)
per= y(35);
cry= y(36);
rev= y(37);
ror= y(38);
bmal= y(39);
Prot_per= y(40);
Prot_cry= y(41);
Prot_rev= y(42);
Prot_ror= y(43);
Prot_bmal= y(44);
PC= y(45);
CB= y(46);
nampt= y(47);
Prot_nampt= y(48);
NAD= y(49);
dbp= y(50);
bmal_prx = y(51);

%% Specify parameter perturbations here
t1=0;
te=5000;
st=0;    %%*((t>t1) & (t<te));
jt=r;    %%*((t>t1) & (t<te));
kt=0;    %%*r*((t>t1) & (t<te));
n=1;     %%*(1+jt);
nj=n;    %%1*(1+jt);
Ki1=Ki;  %%1*(1+kt);


%% Specify LPS dose here
lp=p;
It=lp  ;  %%*10e7*((t>712) & (t<712.1));

%% uncomment this to disconnect HPA axis and inflammation to perform IC50 test
%s0=0;


%%Specifiy dexamethasone dose here
dx=1;  %0
dex = 0 ;
%d=0;

ll=inputF;
if ll == "no"
    for cc = bcz
        if  (t >= (cc+0.1)+d) && (t< (cc+0.26)+d)
            dex=dx*10 ; %Comment this to run IC50 test
        end 
    end 
end

Dex=dex ;

 %Dex=dex*((t>696.00+d) & (t<696.25+d));


%%%%%%%%%%%%%%%% HPA-Axis%%%%%%%%%%%%%
K_strs=10;
Ki_GR1=1.2*Ki ;
n1=1*n;
n2=2;
k_crh=0.096;         
V_acth1=1.4634;%0.9756*1.5;
k_acth=0.4312;%1.15*0.1*2.5*1.5;          
k_cort=0.99;%0.44*0.9*2*1.25;     
V_cort=5.67;%0.0945*2*30;
k_dex1=1.95;
V_crh=0.0667*(1/R);
V_acth2=11.2;
Km_tnf3=40*R;
V_starp2=15;
k_starp=0.45;%0.3*1.5;
V_starp1=0.0225;%0.015*1.5;
Kil6=100;
k_dex2=0.25;
V_dex2=1.15;
Ki_GR=50;
V_GRmrna=3.4;
km_GR1=25*Ki;
kon=0.00329;
k_GRprot=0.0572 ;
krt=0.63 ;
kre=0.57;
k_GRmrna=0.1124;
V_GRprot=1.2;
fr=0.49;
tc=0.1401;
ki_tnf=1e3;
n3=2*n;
n4=2;
fac=16.65;%0.5*33.3;

%parameter values from @ Mavroudis et.al. 2012
kc = 0.009 ; %1/hr    
kac = 5     ; %coupling parameter %1
%IC50rm = 26.2 ;% synthesis drop to half

%%Circadian drive from SCN clock
omega=2*pi/24;
cir=2*(1+ cos(omega*(t)));

%%Circadian cdrive from adrenal peripheral clock
cir2=(4-cir)*(Ki_GR/(Ki_GR+GR)); %(4-cir)*(Ki_GR/(Ki_GR+GR))

%% Corticotrophic releasing hormone (CRH)
dy(1)=(K_strs*(1+st)*(Ki_GR1^n1/(Ki_GR1^n1+GR^n1))*(cir)*(1+V_crh*TNF)-k_crh*CRH);

%% Adrenocorticotropic hormone (ACTH)
dy(2)= (V_acth1*CRH*(Ki_GR1^n1/(Ki_GR1^n1+GR^n1))*(1+V_acth2*((TNF)^n2/(Km_tnf3^n2+ TNF^n2)))-k_acth*ACTH) ;

%% StAR proetin
dy(3)=s0*(V_starp1*(ACTH*cir2)*(1+V_starp2*(TNF^n2/(Km_tnf3^n2+TNF^n2))*(Kil6/(Kil6+IL6)))  -k_starp*StARp);

%% Cortisol
dy(4)= (V_cort*StARp-k_cort*CORT);

%% Dexamethasone kinetics
dy(5)= Dex-k_dex1*Dex1;
dy(6)= V_dex2*Dex1-k_dex2*Dex2;

%% Peripheral cortisol
dy(7)=(1/(tc))*(CORT+(fac*Dex2)-CORTp);

%% GR mRNA
dy(8) = V_GRmrna*(1 - ((GR/(kac*CB))^n3/((km_GR1)^n3+(GR/(kac*CB))^n3)))- k_GRmrna*GR_mrna; 
%dy(8) = V_GRmrna*(1 - (GR^n3/((km_GR1)^n3+GR^n3)))- k_GRmrna*GR_mrna; 

%% GR protein
dy(9) = V_GRprot*GR_mrna + fr*kre*GR- kon*(CORTp)*GR_prot - k_GRprot*GR_prot;

%% GR-cortisol complex in cytosol
dy(10)= kon*(CORTp)*GR_prot- krt*GR_cyt*(ki_tnf^n4/(TNF^n4+ki_tnf^n4)); 

%% Nuclear GR
dy(11)= krt*GR_cyt*(ki_tnf^n4/(TNF^n4+ki_tnf^n4)) -kre*GR;


%%%%% Inflammatory pathway%%%%%%%%%

%%%%%%%Validation parametset 2 for 0.4 ng/kg @ 2 PM%%%%%%%%%%%%%

k_lps=2.7e-5;%1.35*1e-7*100*2;%0
V_phg1= 4.9956*1e7;
V_phg2=12.949;
Km_tnf1=1693.9509;
Ki_tgf1=0.00721;
Ki_IL10=7.384;%147.68*0.05;
k_phg=1.439;
V_tgf1=0.15625*1e-8;
k_tgf=0.0635;%0.03177*2;

V_tnf1=25.5194;
Km_phg1=412500;%0.075*550*1e4;
Ki_tgf2=0.143;%0.15893*0.9;
V_tnf2=106542;% 3.0*3.5514*1e4;
Km_tnf2=123.96;%0.08*1.5495*1e3;
k_tnf=1.25;
Km_il61=80;

Km_phg2=161012;% 8.0506*1e7*0.2*1e-2;
k_il10=1.6;%200*0.008;
V_il102=2.1938e3;%43875*0.05;
Km_tgf3=0.76;% 0.38*2;
V_il101=1.3374e3;%2.67480*1e6*50.0e-5;
Ki_il102=23.636;%1.1818*20*1;
V_il62=5.0e5;
Km_il6tnf=339.164;%4.8452*70;
Km_phg3=11e6;%110*10e4;
V_il61=0.55e5;
k_il6=1.625;%0.5*3.25;
tp=1.5;
tnf_b=1.25;
il6_b=1.5;
V_tgf2=0.5;
Km_GR2=500*Ki1;
n5=1*nj;
n6=4;
n7=2;
n8=2;
n9=6;
n10=2;
cir3=2*(1+ cos(omega*(t)+0.5)) ; % rhythym in peripheral phagocytic cells

%% LPS
dy(12)= 1e-7*(1+It)-k_lps*LPS*Phg;

%% Phagocytes

dy(13)= (V_phg1*(((1*cir3)+ (V_phg2* TNF/(Km_tnf1+TNF)))* (Ki_tgf1/(Ki_tgf1+ TGF))* (Ki_IL10/(Ki_IL10+IL10)))*LPS- k_phg*Phg);

dy(14)=(1.0/(tp))*(Phg-Phg1);

%% Transformig grwth factor (TGF)
dy(15)=  (V_tgf1* Phg + (V_tgf2*GR^n5/(Km_GR2^n5+GR^n5))- k_tgf*TGF);

%% Tumour necrosis factor (TNF)
dy(16)= tnf_b+(Phg/(Km_phg1+Phg))*(V_tnf1+V_tnf2*(TNF/(Km_tnf2+TNF)))*(Ki_tgf2^n6/(Ki_tgf2^n6+TGF^n6))*(1-(IL6^n7/(Km_il61^n7+IL6^n7)))-k_tnf*(TNF);

%% Interleukin 10 (IL10)
dy(17)=((V_il101*Phg1^n8/(Phg1^n8+Km_phg2^n8))+ (V_il102* TGF^n9/(Km_tgf3^n9+ TGF^n9)) - k_il10*IL10);

%% Interleukin 6 (IL6)
dy(18)=il6_b+(Phg1/(Km_phg3+Phg1))*(V_il61+V_il62*(TNF+IL6)^n10/(Km_il6tnf^n10+(TNF+IL6)^n10))*(Ki_il102/(Ki_il102+(IL10)))-k_il6*IL6;



%--------------Meal simulation-------------------

conf=5.31915;%--conversion factor from mg/kg/min to mg/l/min for glucose appearance
%ie-by using distribution volm of glucose=1.88 dl/Kg
Kmaxi=0.0558;%/min
Kmini=0.0080;%/min
Kg_abs=0.057;%%/min------Glu absorption rate connst.
Kp_abs=0.01;%/min---------Protein absn. rate
Kft_abs=0.015;
Kf_abs=0.015;%/min------------Fat absn.rate
Kgri=0.0558;%/min
F=0.90;%dimensionless
%a=0.00013;%/mg
b=0.7;%dim
%c=0.00236;%/mg
d=0.010;%dim

CHO=0;
Pro=0;
Fat=0;

%BW=75;%Kg
%KexA=0.0005;%/min 
%KexB=10;%mmol/l
%b_Qst =0;
%b_Qst1=0;
%b_Qst2=0;
%b_Qgut=0;
%b_Glu_gut=0;
%b_Pro_gut=0;
%b_Fat_gut=0;


%------------Glucose/Protein/Fat Rate of appearance-------

TF=1;
P1=d*0;  
Qst= Qst1+Qst2;
alpha=(2.5./(TF*(1-b)));
  
Kempt= Kmini+(((Kmaxi-Kmini)/2)*((tanh(alpha*(Qst-(b*TF))))+1)); 

%Qst1
dy(19)=-(Kgri*Qst1)+TF*P1;

%Qst2
dy(20)=-(Kempt*Qst2)+(Kgri*Qst1);

%Glu_gut
glf=(CHO/TF);
dy(21)=(Qst2*glf*Kempt)-(Kg_abs*1*Glu_gut);

%Rate Glucose appearance
Ra_Glu=((F*Kg_abs*Glu_gut));%-----converted from mg/min to mmol/min

%Amount of Protein in Gut-Pro_gut
prf=(Pro/TF);
dy(22)=(Qst2*prf*Kempt)-(Kp_abs*Pro_gut);

Kpt_abs=0.01;
dy(23)=(Kp_abs*Pro_gut)-(Kpt_abs*Pro_Int);

%Amount of Fat in Gut Fat_gut
ftf=(Fat/TF);
dy(24)=(Qst2*ftf*Kempt)-(Kft_abs*Fat_gut);

%Rate of appearance of protein-
Ra_Pro=((F*Kpt_abs*Pro_Int)/130);

%Amount of Fat in Gut Fat_gut
ftf=(Fat/TF);

dy(25)=(Kft_abs*Fat_gut)-(Kf_abs*Fat_Int);

%Rate of appearance of Fat-
Ra_Fat=((F*Kf_abs*Fat_Int)/290);%-----converted from mg/min to mmol/min


%---------------------Insulin subsystems---------------
%---Insulin kinetics Parameters--------
Vi=0.05 ;%L/kg
m1=0.190 ;%/min
m2=0.484 ;%/min
m4=0.194 ;%/min
m5=0.0304*Vi ;
m6=0.6471 ;% dimensionless
gamma=0.5;%min

%Insulin secretion
K_AA=0.5;
K_FF=1.8;
na=5.8;
nf=4.8;
n=4.65;
K_Glu=8.9;
V_Glu=48e-12;
V_Ala=17e-12;
V_FFA=21e-12;
Ca_Glu=((Ca_Glup/1.88)*10/180);
Ca_Ala=(Ca_Alap);
Ca_FFA=(Ca_FFAp);
FRnpr= 9.2550;
GCR_Ntv=1.1*(1-0.5* (FRnpr/(20+FRnpr))); % FRnpr = GR-cortisol complex

Ins_sec=V_Glu*(Ca_Glu^n/(Ca_Glu^n+K_Glu^n))...
       +(V_Ala.*((Ca_Ala^na)./(((Ca_Ala^na)+K_AA.^na))))...
       +(V_FFA.*(Ca_FFA.^nf)./((Ca_FFA.^nf)+K_FF.^nf));

ISR=(Ins_sec/1e-12)*GCR_Ntv ;

dy(26)=-(gamma*InsV)+ISR;  % insuline concentration in the portal vien
InSec=gamma*InsV; % Insuline secretion
 
%HE-Hepatic Extraction
HE=((-m5*InSec)+ m6);
m3=((HE*m1)/(1-HE));

%InsL
dy(27)=(-(m1+m3)*InsL)+(m2*InsP)+InSec; % mass of insuline in liver
INSLIV=(InsL/Vi)*1e-12; % insulin concentration in liver

%InsP
dy(28)=(-(m2+m4)*InsP)+ (m1*InsL); % mass of insuline in plasma
INS=1*(InsP/Vi)*1e-12;% insulin concentration in plasma
ki= 0.0079;% 1/min
InsLd = -ki*(InsL-INS);

%Ala pool (release)
Tb_Ala_L=1.5*0.750 ;
Mb_Ala_L=0.25;
Mb_Ala_Lt=0.5;
Ins_P=(1.0*InsP/Vi);
Insu_Ptv_Eff_AA=(10*(Ins_P^2/(Ins_P^2+150^2)));
Reg_AA_up=0.5*(1+Insu_Ptv_Eff_AA);
Ala_up_L=Reg_AA_up*Tb_Ala_L*((Cb_Ala_L/(Mb_Ala_L+Cb_Ala_L))-(0.2/(Mb_Ala_Lt+0.2)));

AA_up_L=Ala_up_L*1.8912;

%Free Fatty Acids(uptake)
Mb_FFA_L=0.57;
Tb_FFA_L=4.77;
Insu_Ptv_Eff_FFA=(10*(Ins_P^3/(Ins_P^3+100^3)));
Reg_Fat_up=0.5*(1+Insu_Ptv_Eff_FFA);

FFA_up_L=Reg_Fat_up*Tb_FFA_L*((Cb_FFA_L/(Mb_FFA_L+Cb_FFA_L))-(0.5/(Mb_FFA_L+0.5))); % need to be discussed bcoz its calculation depends on Blood FFA balance-Cb_FFA

%------Blood metabolite balance----------
Vi=0.05;
Ql=1.5;
INSLIV=InsL/Vi;
INS_b=180e-12;
Meal_bld_Eff=1+(2.5*((INSLIV^4)/((INSLIV^4)+(3*INS_b)^4)));
Qlt=Ql*Meal_bld_Eff;
Vol_tis_L=1.5;
Vol_bld=1.5;
Vol_bld_L=Vol_bld/Vol_tis_L;

dy(29)= (Qlt*(Ca_Alap-Cb_Ala_L)-AA_up_L)/Vol_bld_L; % cb_Ala

dy(30)=(Qlt*(Ca_FFAp-Cb_FFA_L)-FFA_up_L)/Vol_bld_L;

%%%%%%%%% Blood Glucose %%%%%%%
Uii=1*(Ca_Glup^4/(Ca_Glup^4+2^4));%-mg/kg/min
kx1=0.0005;
if Ca_Glup >339
K_ex=kx1*(Ca_Glup-339);
else
K_ex=0;
end

kgp1=0.065;
kgp2=0.079*1.2;

Vm0=2.50*1;%-----mg/kg/min
Vmx=0.047*0.1;%---mg/kg/min
Km0=225.6;%---- mg/kg
Kmx=0;
Vm_X=Vm0+Vmx*Ins_P;
Km_X=Km0;
Ex_eff =1;


Uid_g=(Vm_X*Ca_Glu_OT^1/(Km_X^1+Ca_Glu_OT^1));
 
Ra_Glu1=(Ra_Glu/75);
kp1=2.7; % mg/kg/min
kp2=0.0021;% 1/min
kp3=0.009;% mg/kg/min per pmol/l
kp4=0.0618;% mg/kg/min per pmol/l

EnGlPr = kp1 - kp2*Ca_Glup - kp3*InsLd - kp4*InsV; %Endogenous glucose production

%Plasma Glucose
dy(31)= EnGlPr ...
    + Ra_Glu1...
    -Uii...
    -K_ex...
    -kgp1*Ca_Glup...
    +kgp2*Ca_Glu_OT; % zero represent no change in the plasma glucose  

%Glucose in Other tissues
dy(32)=(-Uid_g+ kgp1*Ca_Glup -kgp2*Ca_Glu_OT);


kap1=0.15*Ca_Alap;
kap2=0.62*Ex_eff;
Cort = 19.5;

Vma0=0.025;
Vmxa=(10*(Ins_P^2/(Ins_P^2+150^2)))*(7/(7+Cort))*1.5;
Kma0=0.6;
Vm_A=Vma0+Vmxa;
Km_A=Kma0;
Uid_p=(Vm_A*Ca_Alap^3/(Km_A^3+Ca_Alap^3));

%Plasma Amino Acids
dy(33)= (-AA_up_L+Ra_Pro-kap1 +kap2-Uid_p); %0;

kfp1=0.1*Ca_FFAp;
kfp2=0.27;%*Ex_eff*arr(124); % Ex_eff , need to be discuss for its value selection either 1 or (0.33*(1+PKA_Ptv_Ex+Cal_Ptv_Ex))

Vmf0=0.025;
Vmxf=1.0*(10*(Ins_P^3/(Ins_P^3+100^3)))*(5*(Cort^2/(Cort^2+7^2)));
Kmf0=1.7;

Vm_F=Vmf0+Vmxf;
Km_F=Kmf0;

Uid_F=(Vm_F*Ca_FFAp^3/(Km_F^3+Ca_FFAp^3));

%Plasma FF Acids
dy(34)= (-FFA_up_L+Ra_Fat-kfp1 +kfp2-Uid_F); % zero represent no chaneg in plasma ffa 0;%




  
%%parameters for mRNA and protein degradatation rates
dm_bmal = 0.827333442085 ;%Bmal mRNA degradation rate
dm_cry = 0.319848706181  ;%Cry mRNA degradation rate
dm_dbp = 0.384486532062  ;%Dbp mRNA degradation rate
dm_nampt = 0.814311309051;%Nampt mRNA degradation rate
dm_per = 0.323114598647  ;%Per mRNA degradation rate
dm_rev = 4.0479072846    ;%Rev-Erb mRNA degradation rate
dm_ror = 0.246575760727  ;%Ror mRNA degradation rate
dp_bmal = 0.186413466797 ;%BMAL protein degradation rate
dp_cry = 0.599026119971  ;%NAMPT protein degradation rate
dp_nampt = 49.8841023982 ;
dp_per = 10.9446515392   ;%PER protein degradation rate
dp_rev = 0.281564063057  ;%REV-ERB protein degradation rate
dp_ror = 0.0340112196281 ;%ROR degradation rate
d_cb =0.197714012552    ;%CLOCK-BMAL complex degradation rate
d_pc = 0.609290138329    ;%PER-CRY complex degradation rate
    
%%parameters for Complexation kinetic rates
Kass_cb = 0.0162923358655 ;%nmol−1· l · h−1 CLOCK-BMAL association rate
Kass_pc = 12.302005485    ;%nmol−1· l · h−1 PER-CRY association rate
Kdiss_cb = 0.00613502224231 ;%h−1 CLOCK-BMAL dissociation rate
Kdiss_pc = 0.0365867175408  ;%h−1 PER-CRY dissociation rate
    
%%parmeters for Maximal transcription rates
Vmax_bmal = 0.0916862770193 ;%Bmal maximal transcription rate
Vmax_cry = 0.702216664807   ;%Cry maximal transcription rate
Vmax_dbp = 0.0802009609453  ;%Dbp maximal transcription rate
Vmax_nampt = 3.49035201578  ;%Nampt maximal transcription rate
Vmax_per = 0.730201742662   ;%Per maximal transcription rate
Vmax_rev = 1.12297601784    ;%Rev-Erb maximal transcription rate
Vmax_ror = 6.9843472736     ;%Ror maximal transcription rate
    
%%parameters for Activation ratios
fold_bmal = 15.9258093373   ;%activation ratio of Bmal by ROR
fold_cry = 1.1604489571     ;%activation ratio of Cry by CLOCK-BMAL
fold_dbp = 400.0            ;%activation ratio of Dbp by CLOCK-BMAL
fold_nampt = 1.57880681573  ;%activation ratio of Nampt by CLOCK-BMAL
fold_per = 12.977351123     ;%activation ratio of Per by CLOCK-BMAL
fold_rev = 73.2342431701    ;%activation ratio of Rev-Erb by CLOCK-BMAL
fold_ror = 335.923333883    ;%activation ratio of Ror by CLOCK-BMAL
    
%%parameters for 
Ka_bmal_ror = 0.00560038941378 ;%Regulation threshold of Bmal by ROR
Ka_cry_cb = 1.0089387144       ;%Regulation threshold of Cry by CLOCK-BMAL
Ka_dbp_cb = 0.308745016237     ;%Regulation threshold of Dbp by CLOCK-BMAL
Ka_nampt_cb = 3.54586790835    ;%Regulation threshold of Nampt by CLOCK-BMAL
Ka_per_cb = 2.03485134763      ;%Regulation threshold of Per by CLOCK-BMAL
Ka_rev_cb = 0.260846828116     ;%Regulation threshold of Rev-Erb by CLOCK-BMAL
Ka_ror_cb = 0.266407416327     ;%Regulation threshold of Ror by CLOCK-BMAL
Ki_bmal_rev0 = 0.0108449480001 ;%Regulation threshold of Bmal by REV-ERB
Ki_cry_rev0 = 0.248955507809   ;%Regulation threshold of Cry by REV-ERB
Ki_cry_pc = 0.00338463577329   ;%Regulation threshold of Cry by PER-CRY
Ki_dbp_pc = 2.23913672671      ;%Regulation threshold of Dbp by PER-CRY
Ki_nampt_pc = 0.0137106537972  ;%Regulation threshold of Nampt by PER-CRY
Ki_per_pc = 0.273493946059     ;%Regulation threshold of Per by PER-CRY
Ki_rev_pc = 28.5885406354      ;%Regulation threshold of Rev-Erb by PER-CRY
Ki_ror_pc = 0.0072858432208    ;%Regulation threshold of Ror by PER-CRY
    
%%parameters for Hill coefficients
hill11 = 4.32985205032; %Hill coeff., regulation of Bmal by REV-ERB
hill10 = 1.83992599578; %Hill coeff., regulation of Bmal by ROR
hill3 = 9.1109447538  ;  %Hill coeff., regulation of Cry by CLOCK-BMAL
hill4 = 2.43715119318 ; %Hill coeff., regulation of Cry by PER-CRY
hill5 = 4.20952050286 ; %Hill coeff., regulation of Cry by REV-ERB
hill14 = 7.32066818222; %Hill coeff., regulation of Dbp by CLOCK-BMAL
hill15 = 10.4312927466; %Hill coeff., regulation of Dbp by PER-CRY
hill12 = 1.91470474775; %Hill coeff., regulation of Nampt by CLOCK-BMAL
hill13 = 1.34080593157; %Hill coeff., regulation of Nampt by PER-CRY
hill1 = 8.52414053707 ; %Hill coeff., regulation of Per by CLOCK-BMAL
hill2 = 8.53897990872 ; %Hill coeff., regulation of Per by PER-CRY
hill6 = 9.83701536835 ; %Hill coeff., regulation of Rev-Erb by CLOCKBMAL
hill7 = 3.31257899336 ; %Hill coeff., regulation of Rev-Erb by PER-CRY
hill8 = 9.36456505724 ; %Hill coeff., regulation of Ror by CLOCK-BMAL
hill9  = 1.84102439743; %Hill coeff., regulation of Ror by PER-CRY
    
%%parameters for Translation rates
Kp_bmal = 0.628507384997 ;%Bmal translation rate
Kp_cry = 3.76912711677   ;%Cry translation rate
Kp_nampt = 58.9647983852 ;%Nampt translation rate
Kp_per = 13.2909782781   ;%Per translation rate
Kp_rev = 0.0513221194724 ;%Rev-Erb translation rate
Kp_ror = 0.0412765888526 ;%Ror translation rate
    
%%parametrs for Protein stability modulation constants
m_cry_ampk = 0.07940386211  ;%modulation of CRY stability by AMPK
m_nampt_ampk = 0.6352243791 ;%modulation of NAMPT stability by AMPK
m_per_ampk = 0.005243953731 ;%modulation of PER stability by AMPK
m_per_sirt = 0.005452322816 ;%modulation of PER stability by SIRT
    
%%parameters for NAD kinetics, Sirt1 and PGC1a activity
Vsirt = 0.915854846427    ;%N/A Maximum SIRT1 activity
Ksirt = 0.75              ;%nmol · l−1 Value of [NAD] at which SIRT1 activity is half of maximum
d_nad = 378.009130749     ;%h−1 Rate of transformation of NAD into NAM
Knad = 0.321746039086     ;%nmol · l−1 Value of [NAD] at which transformationinto NAM is at half of maximum rate
NAD_basal = 0.9116166306  ;%nmol · l−1 Value of [NAD] below which transformation of NAD into NAM is inactivated
Vnad = 296.3933869        ;%molecule per hour per NAMPT protein Maximum regeneration rate of NAD
NAD_tot = 4.166901679     ;%nmol · l−1 Total concentration of NAD and NAM
Knam = 2.76496            ;%nmol · l−1 Value of NAM at which NAD salvage rate is half of maximum
Vpg = 24.06372088         ;%nmol−1· h Maximum activity of PGC1a
Kpg1 = 0.046630145542     ;%N/A Michaelis-Menten constant for phosphorylation of PGC1a by AMPK
Kpg2 = 12.3526351747      ;%N/A Michaelis-Menten constant for deacetylation of PGC1a by SIRT1
    
%%for Pulse parameters 
tc1 = 4.38900149         ;%h Timing of the first AMPK pulse
tc2 = 15.75              ;%h Timing of the second AMPK pulse
tc3 = 18.875             ;%h Time of maximal nuclear PGC1a abundance
Td1 = 2.25               ;%h Duration of the first AMPK pulse
Td2 = 1.5                ;%h Duration of the first AMPK pulse
Td3 = 15.25              ;%h Duration of the nuclear PGC1a presence
Tr1 = 2.6                ;%h Rise time of the first AMPK pulse
Tr2 = 1.8                ;%h Rise time of the second AMPK pulse
Tr3 = 0.5                ;%h Rise time of nuclear PGC1a
amp1 = 6.0               ;%N/A Amplitude of the first AMPK pulse
amp2 = 0.9778008         ;%N/A Amplitude of the second AMPK pulse
amp3 = 0.803062          ;%nmol−1· l Amplitude of the nuclear PGC1a abundance pulse
    
%%paramaters for Chronotherapy timings 
tc4 = 13.664            ;%Timing of the agonist pulse
Td4 = 2.83718           ;%Duration of the agonist pulse
Tr4 = 1.86794           ;%Rise time of the agonist pulse
amp4 = 0.085852;%0.465852;%Amplitude of the agonist pulse

%%Miscellaneous constants used to describe perturbations
%Campk = 1               ;% SIRT1 KO = 1 ; LKB1 KO = 0.0375 ; HFD = 0.05 ; fasting = 0.05
Csirt = 1               ;% SIRT1 KO = 0 ; LKB1 KO = 1 ; HFD = 1 ; fasting = 1
Campk = ampk_calc(Ca_Glup,Ca_Alap,Ca_FFAp);
offs = 0.02             ;% SIRT1 KO = 0.02 ; LKB1 KO = 0.02 ; HFD = 0.02 ; fasting = 2.6
Cpgc1 = 1               ;      

%%calculation of active ampk, sirt1, pgc1a
Prot_PGC1a=  0.28 * per + 0.267 *cry; %amp3 * P3 ;2.35 * per + 2.67 *cry;%
Act_SIRT = (Csirt*Vsirt*NAD)/(Ksirt+NAD) ;


%Act_AMPK =Campk* (amp1*P1 + amp2* P2) + (1 - Campk) * offs ;
Act_AMPK = Campk * (1.57*bmal_prx^2 +  0.07*per^4);
%Campk * (1.356*bmal_prx^2 +  0.125*per^4);%Campk * (1.55*bmal_prx^2 +  0.065*per^4);%

Act_PGC1a= (Cpgc1 * Vpg * Act_AMPK * Act_SIRT * Prot_PGC1a)/(1+((Act_AMPK/Kpg1)*(1+(Act_SIRT/Kpg2)))) ;
agon_rev= amp4 * per^2;%amp4 * P4 ;
Ki_bmal_rev= (Ki_bmal_rev0/(1+agon_rev)) ;
Ki_cry_rev= (Ki_cry_rev0/(1+agon_rev))   ;

cir4= (1+ cos(omega*(t)+7.75));          % 1;

%%for per mRNA estimation 1.38* 
dy(35)= (1.365*Vmax_per*((1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))^hill1)/(1+(CB/(Ka_per_cb*(1+Act_SIRT)))^hill1 * (1+(PC/Ki_per_pc)^hill2))))-dm_per*per + 1.5*kc*(GR/(kac*CB*100)) ;   
%dy(19)= (Vmax_per*((1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))^hill1)/(1+(CB/(Ka_per_cb*(1+Act_SIRT)))^hill1 * (1+(PC/Ki_per_pc)^hill2))))-dm_per*per;  

%%for cry mRNA estimation
dy(36)= (Vmax_cry*((1+fold_cry*(CB/(Ka_cry_cb*(1+Act_SIRT)))^hill3)/(1+(CB/(Ka_cry_cb*(1+Act_SIRT)))^hill3 * (1+(PC/Ki_cry_pc)^hill4) * (1+(Prot_rev/Ki_cry_rev)^hill5))))-dm_cry*cry ;
    
%%for rev mRNA estimation
dy(37)= (Vmax_rev*((1+fold_rev*(CB/(Ka_rev_cb*(1+Act_SIRT)))^hill6)/(1+(CB/(Ka_rev_cb*(1+Act_SIRT)))^hill6 * (1+(PC/Ki_rev_pc)^hill7))))-dm_rev*rev ;
    
%For ror mRNA estimation
dy(38)= (Vmax_ror*((1+fold_ror*(CB/(Ka_ror_cb*(1+Act_SIRT)))^hill8)/(1+(CB/(Ka_ror_cb*(1+Act_SIRT)))^hill8 * (1+(PC/Ki_ror_pc)^hill9))))-dm_ror*ror ;

%%for bmal mRNA estimation 0.41*
dy(39)= (0.44*Vmax_bmal*((1+fold_bmal*(1+Act_PGC1a)*(Prot_ror/Ka_bmal_ror)^hill10)/(1+ (Prot_rev/Ki_bmal_rev)^hill11 + (Prot_ror/Ka_bmal_ror)^hill10)))*cir4-dm_bmal*bmal ;
    
%%for Per protein estimation
dy(40)= Kp_per*per -((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - dp_per*(1+ (m_per_sirt * Act_SIRT) + m_per_ampk * Act_AMPK) * Prot_per ;

%%for Cry protein estimation
dy(41)= Kp_cry*cry -((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - dp_cry*(1+ m_cry_ampk * Act_AMPK) * Prot_cry ;
 
%%for Rev protein estimation
dy(42)= Kp_rev*rev - dp_rev*Prot_rev ; 

%%for Rev protein estimation
dy(43)= Kp_ror*ror - dp_ror*Prot_ror ;
    
%%for bmal protein estimation
dy(44)= Kp_bmal*bmal - ((Kass_cb*Prot_bmal)-Kdiss_cb*CB) - dp_bmal*Prot_bmal ;

%%for Per-Cry complex protein estimation
dy(45)= ((Kass_pc*Prot_cry*Prot_per) - Kdiss_pc*PC) - d_pc* PC ;

%%for Clock-Bmal1 complex protein estimation
dy(46)= (Kass_cb*Prot_bmal - Kdiss_cb*CB) - d_cb* CB ;

%%for nampt mRNA estimation
dy(47)= (Vmax_nampt*((1+fold_nampt*(CB/(Ka_nampt_cb*(1+Act_SIRT)))^hill12)/(1+(CB/(Ka_nampt_cb*(1+Act_SIRT)))^hill12 * (1+(PC/Ki_nampt_pc)^hill13))))-dm_nampt*nampt ;
    
%%for NAMPT Protein estimation
dy(48)= (Kp_nampt*nampt)-((dp_nampt* Prot_nampt)/(1+m_nampt_ampk*Act_AMPK)) ;
    
%%for NAD level estimation
dy(49)= (Vnad*Prot_nampt*((NAD_tot-NAD)/(Knam+NAD_tot-NAD))) - (d_nad*(NAD-NAD_basal)/(Knad+NAD-NAD_basal)) ;

%%for dbp mRNA estimation
dy(50)= (Vmax_dbp*((1+fold_dbp*(CB/(Ka_dbp_cb*(1+Act_SIRT)))^hill14)/(1+(CB/(Ka_dbp_cb*(1+Act_SIRT)))^hill14 * (1+(PC/Ki_dbp_pc)^hill15))))-dm_dbp*dbp ;

dy(51) = 0.55 * bmal- 0.4675*bmal_prx;    
    
df=dy';

end 
