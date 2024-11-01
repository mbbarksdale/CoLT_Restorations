function [dm_equil] = accretion(distance,Ci,amp,numiterations,tr,dt,P,msl,Dmax,Dmin,BMax,ws,timestep,index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use a generic intertidal slope to determine accretion rates accross marsh
% ff=0:.01:1; % flood fraction
C = Ci;
elevation = amp-Dmax:.01:amp;
L = numel(elevation);

for i = 2:numiterations
    depth(i,1:L)=0.5*tr*sin(2*pi*(i*dt/P))+(msl-elevation(1:L));
    ind=find(depth(i,1:L)>0);  % finds indices of points where depth >0
    time_submerged(i,ind)=dt;
    ind2=find(depth(i,1:L)<0);
    time_submerged(i,ind2)=0;
    clear ind
    clear ind2
end

floodfraction=sum(time_submerged,1)/P;

dm = msl + tr/2 - elevation; %gives the depth of the marsh surface below HWL at any given point
for ii = 1:L
    if dm(ii) > Dmax
        bgp(ii) = 0;
    elseif dm < Dmin
        bgp(ii) = 0.0000001; %forest
    else
        bgp(ii)=BMax*(dm(ii)-Dmax)*(dm(ii)-Dmin)/(.25.*(-Dmin-Dmax)*(Dmax-3*Dmin)); %Just refer to biomass here as "bgp" to keep syntax the same, for now.
    end
    bgs(ii) = 0;
end

%vegtype(yr,x)=  % 0=water, 1=scirpus, 2=patens, 3=forest
ind3=find(bgp+bgs==0); vegtype(ind3)=0; bgb(ind3)=0; %water
ind4=find(bgs>bgp);  vegtype(ind4)=1;   bgb(ind4)=bgs(ind4);  %scirpus
ind5=find(bgp>bgs); vegtype(ind5)=2; bgb(ind5)=bgp(ind5);   %patens
mht=msl+(tr/2);

%Now determine mineral deposition
for xx = 2:L
    if(bgb(xx)>0)
        distance(xx)=distance(xx-1)+1;
    else
        distance(xx)=0;
    end
end

for i = 2:numiterations
    ind=find(depth(i,1:L)>0);  % finds indices of points where depth >0
    mineralcycle(i,ind)=C*ws*dt; % no depletion for now
    ind2=find(depth(i,1:L)<0);
    mineralcycle(i,ind2)=0; % no depletion for now
    clear ind
    clear ind2
end

mineral(1:L)=(sum(mineralcycle(:,1:L),1))*timestep; %mineral sedimentation (kg) in an entire year. Summing columns
mineral(1:L)=mineral(1:L).*1000; %Convert from kg to g.

kr=.20; % lignnin content of scirpus, from Ball 1997 
organic(1:L)=kr*bgb*3; %multipled by 3 just to get productivity up reasonably high.
loi(1:L)=organic(1:L)./(mineral(1:L)+organic(1:L))*100;
a=2.67; b=-.7904; %density(yr,x)=a*loi(yr,x).^b; %bulk density, g/cm3, powerlaw fit to Neubauer data
A=loi(1:L); A(A<1)=1;
density(1:L)=(a*A.^b); 
density(1:L)=density(1:L)*100*100*100;
accretion(1:L)=(mineral(1:L)+organic(1:L))./density(1:L); %accretion rate, meters per year

plot(dm,accretion*1000,'Color',[index/131 index/131 index/131])
    
for RSLR = 1:15
    ind = find(accretion*1000>=RSLR,1,'last');
    if isempty(ind) == 0
        dm_equil(RSLR) = amp - elevation(ind);
    else
        dm_equil(RSLR) = Dmax/2;
    end
end
