function [dX,r]=funBAY(X,PAR)

global C_e_ODE
global Fc_ODE
global rhos
global P
global wind
global amp

%Determines change in bay depth and width by solving mass balance between fluxes of sediment into and out of the bay from 
%marsh edge erosion, tidal exchange with the outside sediment source, and sediment deposited onto the marsh surface

%%%%%%%%%%
%defines input parameters for funBAY according to familiar names
ws=PAR(1);
tcr=PAR(2);
Co = PAR(3);
Ba=PAR(4);
Be=PAR(5);
RSLR=PAR(6);
Fm2 = PAR(7);
lamda=PAR(8);
dist=PAR(9);
dmo = PAR(10);
rhob = PAR(11);
rhom = PAR(12);

%%%%%%%%%%%% dynamic variable
fetch=X(1); %mudflat width
df=X(2); %mudflat depth

%%%%%%%%%%%%%%%%%%%%%%%%%
fac=min(1,df/(2*amp));%Proportion of tide that the bay is flooded, dimensionless
Df=(df+(df-fac*2*amp))/2;%[m] average bay depth over tidal cycle
dm = dmo;%[m] marsh edge depth

[tw,Um,Hs,kk] =  wavetau(fetch,wind,Df); %Calculates wave bed shear stress [Pa]

tau=max((tw-tcr)/tcr,0)*lamda;%Excess shear stress, dimensionless
Cr=rhos*tau/(1+tau); %Reference suspended sediment concentration in the basin [kg/m3]

hb=dm+(df-dm)*(1-exp(-dist*0.1/df)); %[m] scarp height at a fixed distance from the marsh according to semi-empirical shoaling profile
W = waveTRNS(amp,wind,fetch,hb); %[W] Wave power density at the marsh boundary

E=(Be*W/(hb-dm)-Ba*Cr*ws/rhom); %(m2/s) Net flux of sediment eroded from/deposited to the marsh edge 

Fc=(Cr-Co)*(fac*2*amp)/P/rhob; %(m2/s) Net flux of sediment lost or gained through tidal exchange with external sediment supply/sink

Fc_ODE(numel(Fc_ODE)+1)=Fc*rhob*fetch; %Save Fc as a mass flux (kg/s) for each iteration of the ODE
C_e_ODE(numel(C_e_ODE)+1) = Cr; %Save C_e (susp. sed. con. at marsh edge, kg/m3) for each iteration of the ODE to use in marsh model
%%%%%%%%%%%%%%%%%%%%%
dX = zeros(2,1); %creates empty matrix for storing outputs from funBAY
dX(1)=E; %(m2/s, or m/s if integrated over 1m transect width) Change in bay width due to erosion
dX(2)=-E*(df-dm)/fetch +Fm2/fetch/rhos + Fc + RSLR; %(m/s) Change in bay depth due to mass balance between fluxes into and out of bay