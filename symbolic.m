%% symbolic matlab

%TOCHECK
% k

%TODO
%eq stato
%lie derivates
%

%% differenziale
syms vby vbx dphi; %stati
syms deltaFL deltaFR; %ingressi

%distanza dal baricentro delle ruote a destra(r) e sinistra(l)
d_l=0.4;
d_r=0.4;

%distanza dal baricentro delle ruote anteriori(f) e posteriori(r)
lf=1;
lr=1;

re=0.25;%raggio ruota

L=lf+lr; %lunchezza viecolo

d=L/tan(deltaFL); %raggio CIR

w_cir=vbx/d;%velocità angolare CIR


%velocità relativa per ogni ruota
v=[sqrt((d-d_l)^2+L^2),...
    sqrt((d+d_r)^2+L^2),...
    d-d_l,...
    d+d_r]...
    *w_cir;

%velocità angolare per ogni ruota
w=v*re;

w_fl=w(1);
w_fr=w(2);
w_rl=w(3);
w_rr=w(4);

%% alphaxy
%side-slip angle approssimati

alphaFL=(vby+dphi*lf)/vbx-deltaFL;
alphaFR=(vby+dphi*lf)/vbx-deltaFR;
alphaRL=(vby-dphi*lr)/vbx;
alphaRR=(vby-dphi*lr)/vbx;


%% k
%slip rate
syms w_FL w_FR w_RL w_RR; %misura velocità ruota


vb=sqrt(vbx^2+vby^2);%velocità complessiva veicolo

%FL
vxFL=vb*cos(alphaFL);%TOCHECK
kFL=-(vxFL-w_FL*re)/vxFL;

%FR
vxFR=vb*cos(alphaFR);
kFR=-(vxFR-w_FR*re)/vxFR;

%RL
vxRL=vb*cos(alphaRL);
kRL=-(vxRL-w_RL*re)/vxRL;

%RR
vxRR=vb*cos(alphaRR);
kRR=-(vxRR-w_RR*re)/vxRR;



%% Gxa and Gyk
%coefficenti per forze trasmesse

Cxa=0.2;
Bxa=0.2;
Exa=0.2;

Cyk=0.2;
Byk=0.2;
Eyk=0.2;

Shyk=0.1;


%Gxa=cos(Cxa*atan(Bxa*alpha-Exa*(Bxa*alpha-atan(alpha))));
%Gyk=(cos(Cyk*atan(Byk*(k+Shyk))))/(cos(Cyk*atan(Byk*Shyk)));

GxaFL=cos(Cxa*atan(Bxa*alphaFL-Exa*(Bxa*alphaFL-atan(alphaFL))));
GykFL=(cos(Cyk*atan(Byk*(kFL+Shyk))))/(cos(Cyk*atan(Byk*Shyk)));

GxaFR=cos(Cxa*atan(Bxa*alphaFR-Exa*(Bxa*alphaFR-atan(alphaFR))));
GykFR=(cos(Cyk*atan(Byk*(kFR+Shyk))))/(cos(Cyk*atan(Byk*Shyk)));

GxaRL=cos(Cxa*atan(Bxa*alphaRL-Exa*(Bxa*alphaRL-atan(alphaRL))));
GykRL=(cos(Cyk*atan(Byk*(kRL+Shyk))))/(cos(Cyk*atan(Byk*Shyk)));

GxaRR=cos(Cxa*atan(Bxa*alphaRR-Exa*(Bxa*alphaRR-atan(alphaRR))));
GykRR=(cos(Cyk*atan(Byk*(kRR+Shyk))))/(cos(Cyk*atan(Byk*Shyk)));


%% Fxy0
%forze normali date dalle caratteristiche dei asfalto, forza normale e
%k/alpha
syms Fz_FL Fz_FR Fz_RL Fz_RR mu_max;

Cx=0.3;
Bx=0.3;
Ex=0.3;

Cy=0.3;
By=0.3;
Ey=0.3;

%FL
Fx0_FL=Fz_FL*mu_max*sin(Cx*atan(Bx*kFL-Ex*(Bx*kFL-Ex*(Bx*kFL-atan(Bx*kFL)))));
Fy0_FL=Fz_FL*mu_max*sin(Cy*atan(By*alphaFL-Ey*(By*alphaFL-atan(By*alphaFL))));

%FR
Fx0_FR=Fz_FR*mu_max*sin(Cx*atan(Bx*kFR-Ex*(Bx*kFR-Ex*(Bx*kFR-atan(Bx*kFR)))));
Fy0_FR=Fz_FR*mu_max*sin(Cy*atan(By*alphaFR-Ey*(By*alphaFR-atan(By*alphaFR))));

%RL
Fx0_RL=Fz_RL*mu_max*sin(Cx*atan(Bx*kRL-Ex*(Bx*kRL-Ex*(Bx*kRL-atan(Bx*kRL)))));
Fy0_RL=Fz_RL*mu_max*sin(Cy*atan(By*alphaRL-Ey*(By*alphaRL-atan(By*alphaRL))));

%RR
Fx0_RR=Fz_RR*mu_max*sin(Cx*atan(Bx*kRR-Ex*(Bx*kRR-Ex*(Bx*kRR-atan(Bx*kRR)))));
Fy0_RR=Fz_RR*mu_max*sin(Cy*atan(By*alphaRR-Ey*(By*alphaRR-atan(By*alphaRR))));


%% Fx and Fy 
%forze risultanti

%FL
Fx_FL=GxaFL*Fx0_FL;
Fy_FL=GykFL*Fy0_FL;

%FR
Fx_FR=GxaFR*Fx0_FR;
Fy_FR=GykFR*Fy0_FR;

%RL
Fx_RL=GxaRL*Fx0_RL;
Fy_RL=GykRL*Fy0_RL;

%RR
Fx_RR=GxaRR*Fx0_RR;
Fy_RR=GykRR*Fy0_RR;


%% equazioni di stato del veicolo tempo continuo

m=1000;%massa veicolo
F_air=0.2*vb;%attrito aria
Izz=250;%inerzia rispetto asse z
%I ECD

dvbx=dphi*vby+...
    (Fx_FL*cos(deltaFL)...
    +Fx_FR*cos(deltaFR)...
    -Fy_FL*sin(deltaFL)...
    -Fy_FR*sin(deltaFR)...
    +Fx_RL+Fx_RR-F_air)/m;

dvby=dphi*vbx+...
    (Fx_FL*sin(deltaFL)...
    +Fx_FR*sin(deltaFR)...
    +Fy_FL*cos(deltaFL)...
    +Fy_FR*cos(deltaFR)...
    +Fy_RL+Fy_RR)/m;

ddphi=(lf*(Fx_FL*sin(deltaFL)+Fy_FL*cos(deltaFL)+Fx_FR*sin(deltaFR)+Fy_FR*cos(deltaFR))...
    -lr*(Fy_RL+Fy_RR)...
    +d_r*(Fx_FR*cos(deltaFR)-Fy_FR*sin(deltaFR)+Fx_RR)...
    -d_l*(Fx_FL*cos(deltaFL)-Fy_FL*sin(deltaFL)+Fx_RL)...
    )/Izz;

syms M_FL M_FR M_RL M_RR;

J=5; %inerzia ruota

% dinamica ruota
dw_FL=(M_FL-Fx_FL*re)/J;
dw_FR=(M_FR-Fx_FR*re)/J;
dw_RL=(M_RL-Fx_RL*re)/J;
dw_RR=(M_RR-Fx_RR*re)/J;

%% equazioni di stato del veicolo tempo discreto

%applicazione Lie-Taylor serie

T0=0.02; %sampling time 

dvbx_dicreto=dvbx+diff(dvbx,vbx)*T0;
dvby_dicreto=dvby+diff(dvby,vby)*T0;
ddphi_dicreto=ddphi+diff(ddphi,dphi)*T0;

dw_FL_dicreto=dw_FL+diff(dw_FL,w_FL)*T0;
dw_FR_dicreto=dw_FR+diff(dw_FR,w_FR)*T0;
dw_RL_dicreto=dw_RL+diff(dw_RL,w_RL)*T0;
dw_RR_dicreto=dw_RR+diff(dw_RR,w_RR)*T0;


