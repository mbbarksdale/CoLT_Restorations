function [DM DB] = wetland3P(Ri,Ci,wind,bfo)

%WETLAND3P.  A 3-point dynamic model for the moprhological evolution of a backbarrier basin with marshes and mudflats
%Copyright (C) 2014, Giulio Mariotti
%Developer can be contacted by <email> and <paper mail>
%This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.function dX=funMARSH(X,PAR)

clearvars -except Ci Ri wind bfo
% close all;

global i
global term
global SSC
global WPD
global erosion

i = 1;

%just to avoid long time close to dfo~0
OPT=odeset('AbsTol',10^-6,'RelTol',10^-6,'Events',@POOLstopp5);
OPTfzero=optimset('Algorithm','Levenberg-Marquardt','TolFun',10^-28,'TolX',10^-28,'MaxFunEvals',10000);

% Ri = 1;
% Ci = 10;

rhos=1000;rhos=2000;rhom=1000;%bulk densities [kg/m3]
%B=10000;%basin width [m]
B=5000;
P=12.5*3600*1; % tidal period [s]
ws=0.5*10^-3; %settlign velocity[m/s]
%ws=0.05*10^-3; %settlign velocity[m/s]
tcr=0.1; %tau cr mudflat [Pa]
amp=1.4/2; % tidal amplitude [m]
Ba=2; %marsh progradation coeff [-]
Be=0.16/(365*24*3600); %marsh erosion coeff [m/yr/ (W/m)]
lamda=0.0001; % mudflat erodability coeff [-]
%dist=3; %reference dist from marsh bank [m].
dist=10;
Co=Ci/1000; % ref conc. [kg/m3]
RSLR=Ri*(10^-3)./(3600*24*365); %relative sea level rise [m/s]

%marsh parameters
BMax=2.500;% kg/m2
nuGp=.0138;
por=1000/2650;
chiref=0.15;

%time
years=100000;TN=1*years+1;to=linspace(1,3600*24*365*years,TN);

filename = ['Outputs/RSLR_' num2str(Ri) '__Co_' num2str(Ci) '/'];

if exist(filename) == 0
    mkdir(filename)
end

%%%%%%%%%%%%%
%OPTION 1, single plot
%%%%%%%%%%%%%%%%

%initial conditions
dfo=2;%initial mudfflat depth
dmo=0.01; %initial marsh depth
% bfo=B/2; %initial mudflat width %Made an input variable

PAR=[rhos P B ws tcr Co wind Ba Be amp RSLR NaN rhom lamda dist bfo BMax nuGp por chiref];
[t,X]=ode23s(@(t,X) funMARSH(X,PAR),to,[bfo dfo dmo],OPT);bf=X(:,1);df=X(:,2);dm=X(:,3);bm=B-bf;


%plot
% Dmax=.7167*2*amp-.483; 
% bf=bf./10000*100;
% axis([0 B./1000 -7.5 1.5])
% plot(bf,-df,'-b',bf,-dm,'-g');hold on; 
% plot(bf(end),-df(end),'ob',bf(end),-dm(end),'dg',[0 B],0*[1 1],'k',[0 B],-Dmax.*[1 1],'k');hold on; 
% plot(bf(1),-df(1),'xr',bf(1),-dm(1),'xr')
% xlabel('b_f [%]')
% ylabel('d_m, d_f [m]')
% axis([0 B/100 -4.5 0.1])
% set(gcf,'units','Inches','position',[1 5 5.8 4.8])

%DCW- plot the depth and width as they vary through TIME
% t_yr = t./(3600*24*365);
% fh2 = figure;
% [hAx,hLine1,hLine2] = plotyy(t_yr,bf,t_yr,df);
% xlabel(hAx(1),'Time (yr)')
% xlabel(hAx(2),'Time (yr)')
% ylabel(hAx(2),'Mudflat depth (m)')
% ylabel(hAx(1),'Mudflatwidth (m)')
% set(hAx(1),'XLim',[0 5000],'YLim',[4800 7700],'YTick',5000:500:7500)
% set(hAx(2),'XLim',[0 5000],'YLim',[1.75 2.65],'YTick',1.8:.1:2.65)

% finalstep = length(bf);
% set(hAx(1),'XLim',[0 finalstep],'YLim',[min(bf)-100 max(bf)+100],'YTick',round(min(bf)/100)*100:500:max(bf))
% set(hAx(2),'XLim',[0 finalstep],'YLim',[min(df)-.2 max(df)+.2],'YTick',round(min(df)/100)*100:.15:max(df))
% 
% text(finalstep*.2,max(bf)-50,['Mudflat depth @ ' num2str(finalstep) 'yr = ' num2str(df(finalstep)) 'm'])
% text(finalstep*.2,max(bf)-250,['Mudflat width @ ' num2str(finalstep) 'yr = ' num2str(bf(finalstep)*100) 'm'])
% set(gcf,'units','Inches','position',[10 5 5.8 4.8])
% 
% %save the plots
% outputfilename2 = [filename 'Bay depth & width plot'];
% print('-dpng',outputfilename2)
% saveas(fh2,[outputfilename2 '.fig'])
% 
% 
% fh3 = figure;
% plot(t_yr,dm)
% xlabel('Time (yr)')
% ylabel('Marsh depth (m)')
% % xlim([0 200])
% text(finalstep*.2,(.01 + dm(finalstep))/2,['Marsh depth @ ' num2str(finalstep) 'yr = ' num2str(dm(finalstep)) 'm'])
% set(gcf,'units','Inches','position',[13 5 5.8 4.8])
% 
% %save the plots
% outputfilename3 = [filename 'Marsh depth plot'];
% print('-dpng',outputfilename3)
% saveas(fh3,[outputfilename3 '.fig'])

% fh4 = figure;
% plot(1:numel(term(1,:)),term(1,:),'r-')
% hold on
% plot(1:numel(term(1,:)),term(2,:),'g-')
% plot(1:numel(term(1,:)),term(3,:),'b-')
% plot(1:numel(term(1,:)),term(4,:),'k-')
% plot(1:numel(term(1,:)),term(5,:),'k--')
% legend('Sediment redistribution from marsh edge','Flux to marsh platform','Flux to ocean','RSLR','Net Change')
% ylabel('Change in bay depth (m)')
% xlabel('Progression through time')
% set(gcf,'units','Inches','position',[1 1 5.8 4.8])
% outputfilename4 = [filename 'Mass balance'];
% print('-dpng',outputfilename4)
% saveas(fh4,[outputfilename4 '.fig'])
% 
% 
% figure
% set(gcf,'units','Inches','position',[6 1 5.8 4.8])
% plot(1:numel(WPD),WPD)
% ylabel('Wave power density (W)')
% figure
% set(gcf,'units','Inches','position',[12 1 5.8 4.8])
% plot(1:numel(erosion),erosion)
% ylabel('Marsh edge erosion (m)')

%%%%%%%%%%%%%
%OPTION 2, multi plot. Recreate Figure 2 in Mariotti and Carr (2014) WRR
%%%%%%%%%%%%%%%%
% CO=[10 50 100]./1000;
% RSLR=[ 5 1 5].*(10^-3)./(3600*24*365);
% in=[.1 .5 .9]; 
%  
% figure
% for i=1:3;
% for j=1:3;bfo=B.*in(j);
%         
% Dmax=.7167*2*amp-.483;                      
% PAR=[rhos P B ws tcr CO(i) wind Ba Be amp RSLR(i) NaN rhom lamda dist bfo BMax nuGp por chiref];
% deq=fsolve(@(X) eq1(X,PAR),2,OPTfzero);
% if deq>0.5;dfo=deq;else;dfo=0.5;end
% dmo=Dmax/2;
% [t,X]=ode23s(@(t,X) funMARSH(X,PAR),to,[bfo dfo dmo],OPT);bf=X(:,1);df=X(:,2);dm=X(:,3);bm=B-bf; 
% 
% subplot(3,1,i)
% bf=bf./1000;
% axis([0 B./1000 -7.5 1.5])
% hold on; 
% plot(bf,-df,'-b',bf,-dm,'-g');hold on; 
% plot(bf(end),-df(end),'ob',bf(end),-dm(end),'dg',[0 B],0*[1 1],'k',[0 B],-Dmax.*[1 1],'k');hold on; 
% plot(bf(1),-df(1),'xr',bf(1),-dm(1),'xr')
% xlabel('b_f [%]')
% ylabel('d_m, d_f [m]')
% axis([0 B/1000 -4.5 0.1])
if numel(dm) >= 200
    DM = dm(1:200);
    DB = df(1:200);
else
    DM = dm;
    DB = df;
end