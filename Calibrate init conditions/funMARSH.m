%WETLAND3P.  A 3-point dynamic model for the morphological evolution of a backbarrier basin with marshes and mudflats
%Copyright (C) 2014, Giulio Mariotti
%Developer can be contacted by <email> and <paper mail>
%This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.function dX=funMARSH(X,PAR)

function dX=funMARSH(X,PAR)

global i
global term
global SSC
global WPD
global erosion

%%%%%%%%%%
%defines input parameters for funBAY according to familiar names
rhos=PAR(1);
P=PAR(2);
B=PAR(3);
ws=PAR(4);
tcr=PAR(5);
Co=PAR(6);
wind=PAR(7);
Ba=PAR(8);
Be=PAR(9);
amp=PAR(10);
RSLR=PAR(11);
rhom=PAR(13);
lamda=PAR(14);
dist=PAR(15);
BMax=PAR(17);
nuGp=PAR(18);
por=PAR(19);
chiref=PAR(20);

%%%%%%%%%%%% dynamic variable
fetch=X(1); %mudflat width
df=X(2); %mudflat depth
dm=X(3); %marsh depth
bm=B-fetch; %marsh width, computed as difference with the basin width

%%%SALT MARSH GROWTH%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Dmin=0; %[m] minimum depth below MHW at which marsh plants can grow
Dmax=(0.237*(amp*2))-0.092+amp; %maximum depth below MHW at which marsh plants can grow, from McKee and Patrick (1988), flexible for different tidal ranges
AA=.25.*(-Dmin-Dmax)*(Dmax-3*Dmin);
Bpeak=BMax*(dm-Dmax)*(dm-Dmin)/AA;
if (Bpeak<=1e-3);Bpeak=0;end
Bfrac=(Bpeak/BMax);
AMC=(180.)*Bpeak*(nuGp)/(365*24*3600);
Rref=AMC*chiref;
FFm=(1/por)*(Rref/rhom);

%%%%%%%%%%%%%%%%%%%%%%%%%
%average depths
fac=min(1,df/(2*amp));
fac2=min(1,dm/(2*amp));
Df=(df+(df-fac*2*amp))/2;
Dm=(dm+(dm-fac2*2*amp))/2; 

tw=  wavetau(fetch,wind,Df,B); %Calculates wave bed shear stress, as a function of fetch, wind speed, and bay depth. From Mariotti and Fagherazzi (2013)
if Dm >1e-4;tw2=wavetauBmod(fetch,wind,Dm,Bfrac);else;tw2=0;end

tau=max((tw-tcr)/tcr,0)*lamda; %Excess shear stress, dimensionless
tau2=max((tw2-tcr)/tcr,0)*lamda;
Cr=rhos*tau/(1+tau); %Reference suspended sediment concentration in the basin [kg/m3]
Cm=rhos*tau2/(1+tau2);
hb=dm+(df-dm)*(1-exp(-dist*0.1/df)); %[m] scarp height at a fixed distance from the marsh according to semi-empirical shoaling profile
W=waveTRNS(amp,df,Df,dist,wind,fetch,dm,hb); %[W] Wave power density at the marsh boundary

E=(Be*W/(hb-dm)-Ba*Cr*ws/rhom); %(m2/s) Net flux of sediment eroded from/deposited to the marsh edge 
Fm=(Cr-Cm)*min(2*amp,dm)/P/rhom; %flux from tidal flat to marsh
Fc=(Cr-Co)*(fac*2*amp)/P/rhom; %(m2/s) Net flux of sediment lost or gained through tidal exchange with external sediment supply/sink

%%%%%%%%%%%%%%%%%%%%%
dX = zeros(3,1);
dX(1)=E;
dX(2)=-E*(df-dm)/fetch +Fm*bm/fetch +Fc +RSLR;

term(1,i) = -E*(df-dm)/fetch;
term(2,i) = Fm*bm/fetch;
term(3,i) = Fc;
term(4,i) = RSLR;
term(5,i) = dX(2);
SSC(i) = Cr;
WPD(i) = W;
erosion(i) = Be*W/(hb-dm);
i = i + 1;


dX(3)=-Fm-FFm+RSLR;