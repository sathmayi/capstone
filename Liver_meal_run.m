clc
clear all


global  n  Ki dex  R p r s0 p1 Q inputF d bcz fFA fAA Glup

inputF= "yes";
p=0.0;%%%LPS injection
r=1.0;%% Change this to vary GR sensitivity
p1=0;%%%% Change this to 1 to introduce LPS for IC50 test
s0=1;%%% Change this to 0 to disconnect HPA with Inflmmation for IC50 test
%u=1;%0.1% make this 0.1 for simulationg IC50 else 1

%% Varying dexamethasone
p_span1=1;%[0,logspace(-2,2.5,20)];%[0:0.01:0.2];%%[5 10 50 100];%[5 10 20 50];% 1.5 2 2.5 3 3.5 4];10;%
%% Varying sensitivity
p_span2=1;%[0.25:0.25:2];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];
%% Varying inhibitory constant
p_span3=1;%[0.1:0.1:1.5];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];

%% Varying cytokine effect on HPA axis
p_span4=1;%[0.1:0.1:1.5];%[0.25:0.25:2];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];

%% Varying GR negative feedback
p_span5=1;%[0.1:0.1:2];%;%[0.5:0.5:4];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];


Ki=p_span3 ;
n=p_span2;
dex=p_span1;
R=p_span4;
Q=p_span5 ;

    
%% the plasma macronutrient levels in fold

fFA = 3;
fAA = 2;
Glup = 2;

x=[];

tspan_max=900;
tspan=linspace(0, tspan_max, tspan_max*10);
options = odeset('RelTol',1e-10,'AbsTol',1.e-10);

bcz = (98:24:tspan_max);

%% Initial conditions
%%%[CRH,ACTH,stARp,CORT, Dex1,Dex2,Cortp,GR_mrna,GR_prot,GR_cyt,GR,LPS,phg,phg1,TGF,TNF,IL10,IL6,per_mrna,cry_mrna,rev_mrna,ror_mrna,bmal_mrna,Prot_per,  Prot_cry ,      Prot_rev   ,   Prot_ror      ,Prot_bmal,      PC,      CB     ,nampt_mrna,  NAMPT,   NAD+, dbp_mrna];
y0=[17.4896364862443 4.82352867112865 0.794555311165795 10.8943173854454 0 0 11.1479571299593 20.1160322148203 342.786973832252...
    24.8798181706611 31.3302385380763 0.000129677801623061 5897.54196590606 25 0.446309929849742 1.00160367257638 22.8093965060946...
    1.09564141359877 0 0 0 0 0 0 0 7.4204 9.4098 2.6370 0.2310 0.5352 166.6464 124.7303 0.3127 0.5679...
    0.552061407767439 1.84289065819700 0.277823248636647 0.994000000000000 1.04468951888596	0.221321027357146...
    1.94826557797849 0.0818595795754448 1.20633513787866 1.55279954517880	8.17708946924250 0.133674457184321 0.791250387033832...
    0.970771596863448	1.12497083683330 0.139125890918912 1.0447];

[t,y] = ode15s(@liver_core_HPA_inflaWmelRes,tspan,y0,options);


%%HPA axis varaiables
CRH=y(:,1);
ACTH=y(:,2);
StARp=y(:,3);
CORT=y(:,4);
Dex1 =y(:,5);
Dex2 =y(:,6);
CORTp =y(:,7);
%%GR regulation related element
GR_mrna =y(:,8);
GR_prot =y(:,9);
GR_cyt =y(:,10);
GR =y(:,11);
%%Inflammation related element
LPS=y(:,12);
Phg=y(:,13);
Phg1= y(:,14);
TGF=y(:,15);
TNF=y(:,16);
IL10=y(:,17);
IL6=y(:,18);
%%From here meal reasponse varaiables starts
Qst1=y(:,19);
Qst2=y(:,20);
Glu_gut=y(:,21);
Pro_gut=y(:,22);
Pro_Int=y(:,23);
Fat_gut=y(:,24);
Fat_Int=y(:,25);
InsV=y(:,26);
InsL=y(:,27);
InsP=y(:,28);
Cb_Ala_L=y(:,29);
Cb_FFA_L=y(:,30);
Ca_Glup	=y(:,31);
Ca_Glu_OT=y(:,32);
Ca_Alap=y(:,33);
Ca_FFAp=y(:,34);
%%From here liver core model strats
per= y(:,35);
cry= y(:,36);
rev= y(:,37);
ror= y(:,38);
bmal= y(:,39);
Prot_per= y(:,40);
Prot_cry= y(:,41);
Prot_rev= y(:,42);
Prot_ror= y(:,43);
Prot_bmal= y(:,44);
PC= y(:,45);
CB= y(:,46);
nampt= y(:,47);
Prot_nampt= y(:,48);
NAD= y(:,49);
dbp= y(:,50);
bmal_prx = y(:,51);

Kpg2 = 12.3526351747 ;
Vpg = 24.06372088    ;
Kpg1 = 0.046630145542;
Vsirt = 0.915854846427;
Ksirt = 0.75          ;


Glup = Ca_Glup(end);
fAA  = Ca_Alap(end);
fFA  = Ca_FFAp(end);

Campk = ampk_calc(Glup, fFA, fAA);

%ampk_calc(Glup,fAA,fFA);

Act_AMPK1=[];
for i = 1:size(t,1)
    Act_AMPK1(end+1) = Campk * (1.57*bmal_prx(i)^2 +  0.07*per(i)^4);
end

Act_PGC1a=[];
for i = 1:size(t,1)
    Act_AMPK = Campk * (1.57*bmal_prx(i)^2 +  0.07*per(i)^4);%Campk * (1.55*bmal_prx(i)^2 +  0.065*per(i)^4);%Campk * (1.356*bmal_prx(i)^2 +  0.125*per(i)^4);
    Prot_PGC1a= 0.28 * per(i) + 0.267 *cry(i);%1.68 * per(i) + 2.72 *cry(i); %amp3 * P3 ;
    Act_SIRT = (Vsirt*NAD(i))/(Ksirt+NAD(i)) ;
    Act_PGC1a(end+1)= (Vpg * Act_AMPK * Act_SIRT * Prot_PGC1a)/(1+((Act_AMPK/Kpg1)*(1+(Act_SIRT/Kpg2)))) ;
end

Ki_bmal_rev=[];
Ki_cry_rev=[];

Ki_bmal_rev0 = 0.0108449480001 ;%Regulation threshold of Bmal by REV-ERB
Ki_cry_rev0 = 0.248955507809   ;%Regulation threshold of Cry by REV-ERB
amp4 = 0.085852;%0.465852;%Amplitude of the agonist pulse

for i = 1:size(t,1)
    agon_rev= amp4 * per(i)^2;%amp4 * P4 ;
    Ki_bmal_rev(end+1)= (Ki_bmal_rev0/(1+agon_rev)) ;
    Ki_cry_rev(end+1)= (Ki_cry_rev0/(1+agon_rev))   ;
end


%%calculation of active ampk, sirt1, pgc1a

figure(1)
hold on;
plot(t,per, '-r');plot(t,bmal, '-g');plot(t,Act_AMPK1, '-b');hold off 

figure(2)
hold on;
plot(t,per, '-r');plot(t,cry, '-g');plot(t,Act_PGC1a, '-b');hold off;  

figure(3)
hold on;
subplot(1,2,1);
plot(t,CORTp, '-r');
subplot(1,2,2);
plot(t,GR, '-g');hold off;  



%{

%---------------Glucose, Amino Acids, Fatty Acids-----------%


figure(1);
hold on;
plot(t,Ca_Glup, '-r');
hold off;

figure(2);
hold on;
plot(t,Ca_Alap, '-g');
hold off;

figure(3);
hold on;
plot(t,Ca_FFAp, '-y');
hold off;



-----------------------------------------
plotting dynamics over 48 hrs
-----------------------------------------

%--------------- HPA variables ---------------%

%cortisol------------

cortisol = y(:, 7);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) - 656;  %start from 0
cort_48h = cortisol(t_range);

figure;
ax = axes;
plot(ax, t_48h, cort_48h, 'b-', 'LineWidth', 2, 'DisplayName', 'Plasma Cortisol (CORTp)');
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1)
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,25,25];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

xlim(ax, [0, 48]);
xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'Cortisol', 'FontSize', 12);
title(ax, 'Cortisol Dynamics over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);


%gr protein----------------

gr_prot = y(:, 9);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) - 656; 

GRprot_48h = gr_prot(t_range);   % GR protein

figure;
ax = axes;
hold on;


dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1)
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,400,400];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, GRprot_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'GR Protein');

xlim(ax, [0, 48]);
ylim(ax, [0, max(GRprot_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'GR protein', 'FontSize', 12);
title(ax, 'GR Protein over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);


% nuclear gr --------------

grn = y(:, 11);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) - 656; 

GRn_48h = grn(t_range);   % GRn

figure;
ax = axes;
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1)
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,40,40];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, GRn_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'nuclear GR');

xlim(ax, [0, 48]);
ylim(ax, [0, max(GRn_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'nuclear GR', 'FontSize', 12);
title(ax, 'nuclear GR over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);


%-----------------Liver Core Clock variables ----------------%

%per gene ---------------


per = y(:, 35);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) - 656; 

per_48h = per(t_range);  

figure;
ax = axes;
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1)
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,40,40];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, per_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'Per mrna');

xlim(ax, [0, 48]);
ylim(ax, [0, max(per_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'per mrna', 'FontSize', 12);
title(ax, 'per over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);


%cry gene --------------

cry = y(:, 36);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) -656;
cry_48h = cry(t_range);  

figure;
ax = axes;
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1) 
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,40,40];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, cry_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'Cry mrna');

xlim(ax, [0, 48]);
ylim(ax, [0, max(cry_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'Cry mrna', 'FontSize', 12);
title(ax, 'Cry mrna over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);


%bmal1 gene --------------

bmal = y(:, 39);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) -656;
bmal_48h = bmal(t_range);  

figure;
ax = axes;
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1) 
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,40,40];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, bmal_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'Bmal mrna');

xlim(ax, [0, 48]);
ylim(ax, [0, max(bmal_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'Bmal mrna', 'FontSize', 12);
title(ax, 'Bmal mrna over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);



% NAD+ ------------------

nad = y(:, 49);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) -656;
nad_48h = nad(t_range);  

figure;
ax = axes;
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1) 
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,6,6];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, nad_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'nad mrna');

xlim(ax, [0, 48]);
ylim(ax, [0, max(nad_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'NAD+', 'FontSize', 12);
title(ax, 'NAD+ over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);



% nampt mrna ------------

nampt = y(:, 47);

t_range = (tspan >= 656) & (tspan <= 704);
t_48h = tspan(t_range) -656;
nampt_48h = nampt(t_range);  

figure;
ax = axes;
hold on;

dark_periods = [12, 24; 36, 48];
labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};
for i = 1:size(dark_periods, 1) 
    start_t = dark_periods(i, 1);
    end_t   = dark_periods(i, 2);
    x_patch = [start_t, end_t, end_t, start_t];
    y_patch = [0,0,6,6];
    patch(ax, x_patch, y_patch, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
end

plot(ax, t_48h, nampt_48h, 'r-', 'LineWidth', 2, 'DisplayName', 'nampt mrna');

xlim(ax, [0, 48]);
ylim(ax, [0, max(nampt_48h) * 1.1]);

xticks(ax, [0, 12, 24, 36, 48]);
xlabel(ax, 'Time (h)', 'FontSize', 12);
ylabel(ax, 'Nampt mrna', 'FontSize', 12);
title(ax, 'Nampt mrna over 48 hours', 'FontSize', 14);
grid(ax, 'on'); ax.GridAlpha = 0.3;
legend(ax, 'show', 'FontSize', 11);


% col_idx, var_name, y_patch_max, color
plot_48h(tspan, y, 7,  'Plasma Cortisol (CORTp)',  25,  'b-');
plot_48h(tspan, y, 9,  'GR Protein',               400, 'r-');
plot_48h(tspan, y, 11, 'Nuclear GR',               40,  'r-');
plot_48h(tspan, y, 35, 'Per mRNA',                 40,  'r-');
plot_48h(tspan, y, 36, 'Cry mRNA',                 40,  'r-');
plot_48h(tspan, y, 39, 'Bmal mRNA',                40,  'r-');
plot_48h(tspan, y, 49, 'NAD+',                     6,   'r-');
plot_48h(tspan, y, 47, 'Nampt mRNA',               6,   'r-');

%}
