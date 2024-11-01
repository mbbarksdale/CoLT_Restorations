function [marshelevation,organic_autoch,organic_alloch,mineral,Fm_min,Fm_org,bgb,accretion,agb] = evolvemarsh(marshelevation,msl,C_e,OCb,sill,SSCreduction,yr,silldelay)

%Calculates biomass and  mineral and organic deposition, for each cell in the marsh as a function of flooding frequency. 
%And calculates the total flux of sediment onto the marsh from the bay.

global tr
global numiterations
global P
global dt
global ws
global timestep
global BMax
global Dmin
global Dmax
global rhoo
global rhos
global SSCreduction
global sill

L = length(marshelevation); %calculate length of the marsh [m]

%Reset variables to zero for the whole marsh domain 
time_submerged(1:L) = 0;
floodfraction(1:L) = 0;
sedimentcycle(1:L) = 0;
depth = 0;
vegtype(1:L) = 0;
bgb(1:L) = 0;
d = 0;

%Loop through a tidal cycle to determine flooding duration for each point in the marsh
for i = 2:numiterations
    depth(i,1:L)=0.5*tr*sin(2*pi*(i*dt/P))+(msl-marshelevation(1:L));%Calculate the depth [m] at each position in the marsh
    ind=find(depth(i,1:L)>0);  % finds indices of points where depth >0
    time_submerged(i,ind)=dt;%These cells are submerged for this portion of the tide cycle
    ind2=find(depth(i,1:L)<0);% finds indices of points where depth <0
    time_submerged(i,ind2)=0;%These cells are not submerged for this portion of the tide cycle
    clear ind
    clear ind2
end

% Calculate belowground productivity. 

%%%%%%%% Biomass curves used by Carr model %%%%%%%%%%%%%%%%%%%%%%
%Creates a biomass curve where peak biomass occurs at a depth halfway between the maximum depth for vegetation and the minimum (here, mean high water level)

dm = msl + tr/2 - marshelevation; %[m] gives the depth of the marsh surface below HWL at any given point
bgb = zeros(1,L); %[g] Belowground Biomass
agb = zeros(1,L); %[g] Aboveground Biomass 
organic_autoch=zeros(1,L); %creates empty matrix for storing autochthonous OM values

for ii = 1:L %starts for-loop running through length of the marsh
    if dm(ii) > Dmax %If depth is below vegetation maximum, there is no production
        agb(ii) = 0; %[g] mudflat
        bgb(ii) = 0; %[g] mudflat
        organic_autoch(ii)=0; %[g] no autochthanous material stored in the soil;
    elseif dm(ii) <= Dmin  %If depth is above vegetation minimum, there is very little belowground productivity
        if ii>=5500 %if marsh position is greater than initial forest edge ... [not flexible for different marsh widths!]
            agb(ii)=100; %[g]forest aboveground - low background level
            bgb(ii) = 0.00001; %[g] forest
            organic_autoch(ii)= bgb(ii); %forest organic matter stored in soil %[g] 
        else %if marsh position is not greater than initial forest edge ...
            agb(ii) = 1*(BMax*(dm(ii)-Dmax)*(dm(ii)-Dmin)/(.25.*(-Dmin-Dmax)*(Dmax-3*Dmin))); %aboveground biomass follows biomass curve
            bgb(ii) = 0.1; %[g] "forest" belowground biomass - really high marsh
            organic_autoch(ii)= bgb(ii); %forest organic matter stored in soil %[g] we do not multiply by lignin content (as in earlier version of the model) because we subtract mass due to decomposition in the function 'decompose'
        end
    else %if depth below mean high water is between Dmin and Dmax, calculate according to the relationship in Mariotti and Carr 2014, based on Morris et al 2002
        agb(ii)=max(1*(BMax*(dm(ii)-Dmax)*(dm(ii)-Dmin)/(.25.*(-Dmin-Dmax)*(Dmax-3*Dmin))),0); %[g] aboveground biomass follows marsh biomass curve
        bgb(ii)=max(BMax*(dm(ii)-Dmax)*(dm(ii)-Dmin)/(.25.*(-Dmin-Dmax)*(Dmax-3*Dmin)),0); %[g] belowground biomass follows marsh biomass curve
        organic_autoch(ii)=bgb(ii); %[g] we do not multiply by lignin content (as in earlier version of the model) because we subtract mass due to decomposition in the function 'decompose'
    end
end

%Marsh width only includes those cells with belowground biomass
marshwidth=length(find(bgb>0)); %[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine suspended sediment concentration at marsh edge with sill function
if sill==0 || yr < silldelay %if sill is turned off in input values OR if there is a delay in sill implementation/installation
    C_e=C_e; %SSC at marsh edge doesn't change
else if sill==1 && yr >= silldelay %if sill is turned on and any delays in sill implementation/installation occurred before the current year... 
    C_e=C_e*(100-SSCreduction)*0.01; % this decreased SSC is based on decrease in water velocity behind the sill, which would increase the falling out of suspended particles
else %if none of these conditions are met...
    error('Error. Must choose either 0 or 1 for Input file column 6.') %display this error message
end
end

%Now determine mineral deposition
distance = 0; %[m] %sets distance marker to 0 at start and will increase during for-loop
for xx = 1:L %for the length of the marsh
    if bgb(xx)>0 %if there is vegetation
        distance=distance+1; %[m] distance from marsh edge increases by 1 m
        C(xx) = C_e*exp(-0.0031.*distance); %[kg/m3] calculate suspended sediment concentration across the marsh (exponential decay function), coefficient "-.0031" is a fitted parameter for realistic marsh topography
    else %if there is no vegetation
        distance=1; %[m] reset distance from marsh edge
        C_e=0.5*C_e; %[kg/m3] Decrease concentration at the new "marsh edge" (i.e. the pond) by 50% with each subsequent pond formation, parameterizing channels importing sediment into ponds
        %C(xx) = C_e*exp(-0.0031.*distance); %[kg/m3] coefficient "-.0031" is a fitted parameter for realistic marsh topography
        C(xx) = C_e; %[kg/m3] suspended sediment concentration in the pond
    end
end

floodfraction=sum(time_submerged,1)/P; %Portion of the tidal cycle that each point is submerged

for i = 2:numiterations
    ind=find(depth(i,1:L)>0);  % finds indices of points where depth >0
    sedimentcycle(i,ind)=C(1,ind)*ws*dt; %[kg] mass of mineral sediment deposited
    ind2=find(depth(i,1:L)<0);% finds indices of points where depth <0
    sedimentcycle(i,ind2)=0; %[kg] no mineral sediment deposited
    clear ind
    clear ind2
end

susp_dep=sum(sedimentcycle(:,1:L),1)*timestep; %[kg/yr] suspended sediment deposition in an entire year, in each cell (timestep = # tidal cycles/yr)
susp_dep=susp_dep.*1000; %[g] Convert from kg to g.

%plot sediment dynamics every year
figure(100)
subplot(3,1,1)
plot(sedimentcycle')
subplot(3,1,2)
plot(C)
hold on
subplot(3,1,3)
plot(susp_dep)
hold on

mineral=susp_dep.*(1-OCb); %[g] Mineral deposition of suspended sediment in a given year is determined by the mineral content of the bay sediment 
organic_alloch=susp_dep.*OCb; %[g] Organic deposition of suspended sediment is determined by the organic content of the bay sediment 

Fm_min = sum(mineral)/1000; %[kg/yr] Flux of mineral sediment from the bay
Fm_org = sum(organic_alloch)/1000; %[kg/yr] Flux of organic sediment from the bay

Fmin = mineral(1); %Fmin is defined as mineral deposition of first marh cell
Fp = 0; %Flux of sediment from marsh pond that forms when accretion in the interior marsh is not sufficient to keep pace with SLR
recharge = 0; % Signifies whether sediment concentration has been recharged by pond formation
for xxx = 2:L %starts for-loop to run through length of marsh domain
    if C(xxx) < C(xxx-1) && recharge == 0 %as long as sediment concentration is decreasing, and hasn't been recharged by pond formation, sediment flux is coming from bay
        Fmin = Fmin + mineral(xxx); %Flux of mineral sediment is current flux plus mineral  
    else %if sed concentration isn't decreasing (which happens when ABG=0), or if recharge = 1, then every position after ponding receives pond sed + min dep
        recharge = 1;
        Fp = Fp + mineral(xxx); %once pond formation results in sediment recharge, sediment flux comes from pond
    end
end

% Calculate thickness of new sediment (mineral+organic) based off LOI and its effect on density
loi=(organic_autoch+organic_alloch)./(mineral+organic_autoch+organic_alloch); %[LOI]
density=1./((loi./rhoo)+((1-loi)./rhos)); %[kg/m3] Bulk density is calculated according to Morris et al. (2016)
density=density.*1000; %[g/m3]
density(isnan(density))=1; %If there is zero deposition, loi calculation will divide by zero and make density nan. Set density equal to one in this case, so that accretion is zero, instead of nan.

accretion(1:L)=(mineral+organic_autoch+organic_alloch)./density(1:L); %[m] accretion in a given year

%loi and density aren't changing through time because all cells are being replaced each time interval.

%Update elevation
marshelevation=marshelevation+accretion; %[m]