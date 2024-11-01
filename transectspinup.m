function transectspinup(filename,runi)

%transectspinup. A dynamic model spinup for the morphological evolution of a backbarrier basin with marshes, mudflats, and an upland slope, and their requisite carbon pools.

close all;clc;clear global

global C_e_ODE %edge erosion model
global Fc_ODE %model calculating flux of sediment into and out of bay
global tr %[m] tidal range
global numiterations %number of increments in tidal cycle
global P %[s] tidal period
global dt %[s] seconds per numiteration
global ws %[m/s] effective settling velocity
global timestep %[tidal cycles per year] number to multiply accretion simulated over a tidal cycle by
global filename %name of file where run data will be stored
global endyear %last year of model run

global Dmin %[m] minimum depth below MHW at which marsh plants can grow
global Dmax %[m] maximum depth below MHW at which marsh plants can grow
global BMax %[g/m^2] maximum plant biomass, at optimal conditions
global amp %[m] tidal amplitude
global wind %[m/s] wind speet
global bfo %[m] initial bay fetch -- If this is changed, dfo and dmo should be recalibrated
global B %[m] length of domain

global fetch %[m] length of fetch/mudflat/bay
global yr %[yr] current year of experiment
global mui %[m] depth below which decomposition goes to zero in the marsh
global mki % Coefficient of decomposition in the marsh
global rhoo %[kg/m3] bulk density of organic matter 
global rhos %[kg/m3] bulk density of mineral matter 
global Coi %[mg/L] external sediment supply

Inputdata=xlsread(['transectspinupfiles\' filename '\Input variables.xlsx']); % Read in input file, which is an excel file located in the "Transect Spinup Files" --> "spinups" folder within the larger folder containing the model

%this section reads in the input files info from excel sheet. Read in the input files info from excel sheet and rename the variable to desired names
if nargin == 1 % Case 1: Running only 1 model run
    RSLRi=Inputdata(1,1); %[mm/yr] relative sea level rise rate at beginning of spinup (i=initial)
    RSLRf=Inputdata(1,2); %[mm/yr] relative sea level rise rate at end of spinup (f=final)
    Coi=Inputdata(1,3); %[mg/L] external sediment supply
    slope=Inputdata(1,4); %[m/m]Upland slope
    dur=Inputdata(1,5); %[yr] Number of years to simulate
else  %Case 2: Running a batch of model runs
    RSLRi=Inputdata(runi,1); %[mm/yr] relative sea level rise rate at beginning of spinup (i=initial)
    RSLRf=Inputdata(runi,2); %[mm/yr] relative sea level rise rate at end of spinup (f=final)
    Coi=Inputdata(runi,3); %[mg/L] external sediment supply
    slope=Inputdata(runi,4); %Upland slope
    dur=Inputdata(runi,5); %[yr] Number of years to simulate
end

%Convert variables to more useable units    
SLRend = RSLRf; %saves another variable as SLRend
RSLR(1:299)=RSLRi; %establishes the first 299 years as the initial SLR
for y=300:400 %begins a loop 
    RSLR(y) = ((SLRend-RSLRi)/100)*(y-300) + RSLRi; %establishes the rate of RSLR during years 300 to 400 according to a linear ramp between the initial and final rates of SLR 
end
RSLR(401:dur)=RSLRf; %establishes the remaining years of the spinup as the final rate of SLR, allowing marsh accretion rates to approach equilibrium
RSLR=RSLR*(10^-3)./(3600*24*365); %relative sea level rise rate [m/s]
SLR=RSLR*(3600*24*365); %Sea level rise for a given yr [m/yr]

Co=Coi/1000; % ref conc. [kg/m3] convert reference concentration to kg/m3 (from mg/L in input sheet)

%create empty vectors for OM storage
organic_dep_autoch=zeros(dur,B);
organic_dep_alloch=zeros(dur,B);
bgb=zeros(dur,B);
agb=zeros(dur,B);

startyear=2; %force starting elevation for year 1, so startyear is the year that marsh spinup beings
endyear = dur; %Last year to simulate 

%Name of the folder where outputs/plots will be saved
outputfilename = ['transectspinupfiles/' filename '/Spinup_CO' num2str(Coi) '_RSLRinitial' num2str(RSLRi) '_rampedRSLRfinal' num2str(RSLRf) '/'];
if exist(outputfilename,'file') ~= 7 %if filenmae does not exist...
    mkdir(outputfilename) %...create a folder to save the output variable to
end

OPT=odeset('AbsTol',10^-6,'RelTol',10^-6,'Events',@POOLstopp5); %set tolerances for ODE solver
OPTfzero=optimset('Algorithm','Levenberg-Marquardt','TolFun',10^-28,'TolX',10^-28,'MaxFunEvals',10000); %set tolerances for ODE solver

%hard-coded parameters
rhos=2000;%bulk density of mineral matter [kg/m3]
rhob=rhos;%bulk density of the bay, which is initially all minerl [kg/m3]
rhoo=85.0;%bulk density of organic matter [kg/m3]
P=12.5*3600*1; % tidal period [s]
ws=0.05*10^-3;%effective settling velocity [m/s]; Marani 2010 uses 0.2 mm s-1, Mudd 2009 uses .037 mm/s ws2=0.5*10^-3;
wsf=0.5*10^-3; %settling velocity on the mud flat[m/s]
tcr=0.15; %tau cr mudflat [Pa] 
wind=6; %reference wind speed [m/s]
amp=0.5; % tidal amplitude [m] 
Ba=0; %no progradation of marsh during spinup
Be=0; %no erosion of marsh during spinup
lamda=0.0001; % mudflat erodability coeff [dimensionless] (Mariotti and Carr, 2014)
dist=10; %reference dist from marsh bank [m]. (Mariotti and Carr, 2014)
bfo=5000; %initial bay fetch -- If this is changed, dfo and dmo should be recalibrated

tr=amp*2; %[m] tidal range
timestep=365*(24/12.5); %[tidal cycles per year] number by which to multiply accretion simulated over a tidal cycle

msl(1)=0; %sets mean sea level at 0 for starting year
for y=2:endyear
    msl(y)=SLR(y)+msl(y-1); %[m] mean sea level, calculates MSL over time given the SLR rate
end

%Minimum and maximum depth below high water that marsh vegetation can grow
% K the standard parabola
BMax=2500;% [g/m2] maximum biomass production for marsh vegetation  
Dmin=0; %[m] minimum depth below MHW at which marsh plants can grow
Dmax=(0.237*tr)-0.092+amp; %maximum depth below MHW at which marsh plants can grow, from McKee and Patrick (1988), flexible for different tidal ranges

%Variables for forest aboveground growth
BMax_forest=5000; %[g/m2] maximum biomass production for forest vegetation  
a=4; %forest growth coefficient 
b=2; %Tree growth rate [1/m]
f0=0.0001; %Background carbon accumulation in forest soils [g/m2]
fwet=5; %Carbon layer from wetted soils [g/m2]
fgrow=2; %Decay constant for soil carbon [1/m]

%Set up domain and time vectors
years=1;TN=1*years+1;to=linspace(1,3600*24*365*years,TN); %initiate time, create vector for time
x_m = ceil(bfo)+1; %first marsh cell
x_m = ceil(bfo)+1; %Marsh edge (first marsh cell)
Marsh_edge = zeros(endyear,1); %create vector to track marsh edge
Marsh_edge(1:startyear-1)=x_m; %designate marsh edge position during the spinup (no erosion during spinup)
Forest_edge=zeros(endyear,1); %create vector to track forest edge

% Calculate inundation time and mineral sedimentation in a single tidal cycle
numiterations=500; %number of increments by which divide tidal cycle in evolvemarsh function
dt=P/numiterations; %[s] seconds per numiteration (here 500 in a single tidal cycle and thus dt=90 seconds per timestep) 

OCb=zeros(1,endyear); %Organic content of the uppermost layer of bay sediment, which determines the OC of suspended material deposited onto the marsh. Initially set to zero.
OCb(1:550)=.05; %organic content of bay for spinup years set to 5%
edge_flood = zeros(1,endyear); %create vector to track marsh drowning
Fm_min = 0; %[kg/yr] sets the flux of mineral sediment from the marsh to the bay to 0
Fm_org = 0; %[kg/yr] sets the flux of organic sediment from the marsh to the bay to 0
Fm_flood = 0; %[kg/yr] sets the flux of sediment from the marsh to the bay due to flooding to 0
Fp_sum = 0; %The amount of sediment taken from ponds to recharge sedimentation to drowning interior marsh
fetch = zeros(1,endyear); %create vector to track fetch throughout simulation
fetch(1) = bfo; %[m] set fetch to initial value for duration of the spinup
x_b=1; %sets initial first bay cell

%Populate first 25 years of deposition with underlying stratigraphy. To change initial marsh morphology or stratigraphy, it must be done in the "Spin Up" file.
mwo=500; %Initial marsh width

% if changing the spinup to have sudden SLR change, need to change several lines of code in buildtransect function!
[B,db,elevation] = buildtransectspinup(amp,wind,bfo,filename,endyear,RSLRi,msl,Coi,slope,mwo); %Build initial transect for spinup, see buildtransectspinup function for more detail
Forest_edge(startyear-1)=find(elevation(startyear-1,:)<(msl(startyear-1)+amp-Dmin),1,'last')+1; % % Calcuate where elevation is right for the forest to start

%K set up vectors for deposition
organic_dep_alloch=zeros(endyear,B);
organic_dep_autoch=zeros(endyear,B);
mineral_dep=zeros(endyear,B);

Bay_depth(1:startyear) = db; %aligns bay depth with db from buildtransectspinup function
dmo = elevation(startyear-1,x_m); %Set marsh edge depth to the elevation of the marsh edge at end of spinup
mortality = zeros(endyear,B); %Create vector for plant mortality
BayExport = zeros(endyear,2); %Create vector for bay export
BayOM = zeros(endyear,1); %Create vector for OM in the bay
BayMM = zeros(endyear,1); %Create vector for minearal matter in the bay
fluxes = zeros(8,endyear); %Create matrix to store fluxes

%Decomposition parameters
mui = 0.4; %[m] Depth below which decomposition goes to zero in the marsh
mki = 0.1; %Coefficient of decomposition in the marsh

bgb_sum=zeros(endyear,1); %[g] set up vector for sum of organic matter deposited across the marsh platform in a given year
Fd=zeros(endyear,1); %[kg] set up vector for flux of organic matter out of the marsh due to decomposition
avg_accretion=zeros(endyear,1); %[m/yr] set up vector for annual acretion rate averaged across the marsh platform

%Load in lookup table with soil organic matter for forest based on age and depth
load('Forest Organic Profile/forestOM.mat') 
load('Forest Organic Profile/forestMIN.mat') 
load('Forest Organic Profile/B_rts.mat') 

forestOM; %[g] Table with forest organic matter profile stored for in 25 depth increments of 2.5cm (rows) for forests of different ages (columns) from 1 to 80 years
forestMIN; %[g] Table with forest mineral matter profile stored for in 25 depth increments of 2.5cm (rows) for forests of different ages (columns) from 1 to 80 years
startforestage=60; %[yr] age of the forest at the start of the simulation
forestage=startforestage; %set forest age to the forest age specified for the beginning of the simulation

%start a time loop, where you run funBAY every ts to determine change in bay depth and width, then run the evolvemarsh to determine marsh elevation change, then go to next time step and iterate with updated values.
for yr = startyear:endyear 
    
    yr % print out year of iteration to see progress while running
    
    %Calculate the density of the marsh edge cell
    boundyr=find(elevation(:,x_m)>elevation(yr-1,1),1); %Find first year where the elevation of the marsh is above the depth of erosion (i.e. bay bottom)
    if boundyr == 1 %If bay bottom is not eroded...
        us=elevation(1,x_m)-elevation(yr-1,1); %[m] Depth of underlying stratigraphy
        usmass=us*rhos; %[kg] Mass of pure mineral sediment underlying marsh at marsh edge
    else
        usmass=0; %[kg] Mass of pure mineral sediment underlying marsh at marsh edge
    end
   
    massm=sum(organic_dep_autoch(boundyr:yr-1,x_m))/1000 +sum(organic_dep_alloch(boundyr:yr-1,x_m))/1000 + sum(mineral_dep(boundyr:yr-1,x_m))/1000 + usmass; %[kg] Mass of sediment to be eroded at the current marsh edge above the depth of erosion
    volm=elevation(yr-1,x_m)-elevation(yr-1,1);  %[m3] Volume of sediment to be eroded at the current marsh edge above the depth of erosion
    
    rhom=massm/volm; %bulk density of marsh edge [kg/m3]
    if rhom>rhos %if bulk density of marsh edge is greater than pure sediment
        rhom=rhos; %rhom should never exceed the density of pure sediment
    end
    if rhom<rhoo %if bulk density of marsh edge is less than density of pure OM 
        rhom=rhoo; %rhom should never be less than the density of pure OM
    end
    rhomt(yr)=rhom; %create vector to track density of eroded material over time (saved at end of simulation as output)
    massmt(yr)=massm; %create vector to track mass of eroded material over time (saved at end of simulation as output)
    Fm=(Fm_min+Fm_org)/(3600*24*365); %[kg/s] Mass flux of both mineral and organic sediment from the bay to the marsh
    
    Erosion=zeros(1,endyear); %create vector to save erosion through time

    PAR=[wsf tcr Co Ba Be RSLR(yr) Fm lamda dist dmo rhob rhom]; % Set up parameters to feed into ODE
    [t,X]=ode23s(@(t,X) funBAY(X,PAR),to,[bfo db],OPT);fetch_ODE=X(:,1); db_ODE=X(:,2); %Solves for change in bay depth and width, solver for still diff equat., low order method, see funbay for details of the calcuation

    Erosion(yr)=0;

    %Remove any NaN values from bay fetch and depth results
    db_ODE(isnan(db_ODE))=[];
    fetch_ODE(isnan(fetch_ODE))=[];
    Fc_ODE(isnan(Fc_ODE))=[];
    
    %update values for this year iteration
    db = db_ODE(end); %Set initial bay depth of the bay to final depth from funBAY
    fetch(yr) = fetch_ODE(end); %Set initial bay width of the bay to final width from funBAY
    baydepth=db; %aligns baydepth variable with db (both used)
    bfo = fetch(yr); %Set initial bay width of the bay to final width from funBAY
    C_e(yr) = C_e_ODE(end);%[kg/m^3] Susp. sed. conc. at marsh edge should be defined as an output from  edge erosion model
    Fc=Fc_ODE(end)*(3600*24*365); %[kg/yr] Annual net flux of sediment out of/into the bay from outside the system %MBB see code above and can delete if not useful! 
    Fc_org=Fc*0.05; %[kg/yr] Annual net flux of organic sediment out of/into the bay from outside the system, based on a 5% organic load in incoming sediment
    Fc_min=Fc*0.95; %[kg/yr] Annual net flux of mineral sediment out of/into the bay from outside the system, based on a 95% mineral load in incoming sediment
    
    [Fe_org,Fe_min] = calcFE(bfo,fetch(yr-1),baydepth,elevation,organic_dep_autoch,organic_dep_alloch,mineral_dep,msl,boundyr); %Calculate the flux of organic and mineral sediment to the bay from erosion of the marsh
    Fe_org=Fe_org/1000; %[kg/yr] Annual net flux of organic sediment to the bay due to erosion
    Fe_min=Fe_min/1000; %[kg/yr] Annual net flux of mineral sediment to the bay due to erosion
    
    Fb_org = Fe_org - Fm_org - Fc_org; %[kg/yr] Net flux of organic sediment into (or out of, if negative) the bay
    
    Fb_min = Fe_min - Fm_min - Fc_min; %[kg/yr] Net flux of mineral sediment into (or out of, if negative) the bay
    
    BayExport(yr,:) = [Fc_org Fc_min]; %[kg/yr] Mass of organic and mineral sediment exported from the bay each year
    BayOM(yr)=Fb_org; %[kg/yr] Mass of organic sediment stored in the bay in each year
    BayMM(yr)=Fb_min; %[kg/yr] Mass of mineral sediment stored in the bay in each year
    
    
    if Fb_org > 0 && Fb_min > 0 %If there is both mineral and organic flux into the bay
        OCb(yr) = Fb_org/(Fb_org+Fb_min); %[%] bay organic content for each year is defined as contribution of organic flux into bay divided by total flux into bay
    elseif Fb_org > 0 %if there is only organic flux into the bay, organic content of the bay bottom surface =100%
        OCb(yr) = 1;
    elseif Fb_min > 0 %if there is only mineral flux into the bay, organic content of the bay bottom surface = 0%
        OCb(yr) = 0;
    else %If there is no mineral or organic deposition, the bay bottom surface has the same organic content as the previous timestep
        OCb(yr) = OCb(yr-1);
    end
    
    if db > Bay_depth(1) %If bay has eroded down to depth below initial bay bottom, there is only mineral sediment remaining
        OCb(yr) = 0.05;
    end
    
    rhob = 1 / ((1-OCb(yr))/rhos + OCb(yr)/rhoo); % [kg/m3] Density of bay sediment
    
    if isreal(bfo) == 0 %if the fetch is not a real number, means that the marsh is completely gone and the simulation ends! 
        disp('Marsh has eroded completely away')
        endyear = yr;
        break
    end
    
    x_m = ceil(bfo)+1; %new first marsh cell (i.e., marsh edge)
    %x_m = find(elevation(yr-1,:)>0,1); %first new marsh cell (i.e., marsh edge)
    x_f = find(elevation(yr-1,:)<(msl(yr)+amp-Dmin),1,'last')+1; %new first forest cell (i.e., forest edge)
       
    tempelevation = elevation(yr-1,x_m:x_f-1); %for ease off coding, create tempelevation variable to use in calcuations each year (previous year's elevation)
    Dcells = Marsh_edge(yr-1)-x_m; %Calculates change in the marsh edge from previous year to current year (x_m) 
    
    if Dcells > 0 %Prograde the marsh, with new marsh cells having the same elevation as the previous marsh edge
        tempelevation(1:Dcells) = elevation(yr-1,Marsh_edge(yr-1)); %Prograde the marsh, with new marsh cells (from marsh edge to previos marsh edge) having the same elevation as the previous marsh edge
    end
    
    msl(yr) = msl(yr-1) + SLR(yr); %updates msl for current year
    elevation(yr,1:x_m-1) = msl(yr)+amp-db; %All bay cells have the same depth, set to the new value, taking into account the MSL

    [tempelevation,temporg_autoch,temporg_alloch,tempmin,Fm_min,Fm_org,tempbgb,accretion,tempagb] = evolvemarshspinup(tempelevation,msl(yr),C_e(yr),OCb(yr)); %allow for deposition in the marsh, see function evolvemarsh for more detail

    elevation(yr,x_m:x_f-1) = tempelevation; %[m] update elevation of the given year from the tempelevation that has gone through evolvemarsh   
    elevation(yr,x_f:B) = elevation(yr-1,x_f:B); %Forest elevation remains unchanged
    mineral_dep(yr,x_m:x_f-1) = tempmin; %[g] mineral sediment deposited in a given year    
     
    %AUTOCHTHONOUS WITH DEPTH
    dddd=elevation(yr,x_m:x_f-1) - elevation(1:yr,x_m:x_f-1); %Calculate depth of sediment pocket below the surface
    [CC]=biodepth(dddd,tempbgb,x_m,x_f-1,organic_dep_autoch,elevation); %biodepth function distributes belowground biomass to depths of up to 0.4 m (refer to function for more details) 
    organic_dep_autoch=organic_dep_autoch+CC; %add in new OM to previous year's OM

    mortality(yr,x_m:x_f-1) = temporg_autoch; %[g] belowground plant material deposited in a given year, for keeping track of without decomposition
    organic_dep_alloch(yr,x_m:x_f-1) = temporg_alloch; %[g] allochthanous organic material deposited in a given year
    bgb_sum(yr) = sum(tempbgb); %[g] belowground biomass deposition summed across the marsh platform. Saved through time without decomposition for analysis
    %Fm_min=Fm_min+Fm_min_prog; %[kg/yr] Add fluxes deposited at marsh edge to fluxes deposited on marsh platform
    %Fm_org=Fm_org+Fm_org_prog; %[kg/yr] Add fluxes deposited at marsh edge to fluxes deposited on marsh platform
    bgb(yr,x_m:x_f-1)=tempbgb; %[g] adds the year's marsh belowground biomass to the matrix in its correct row
    agb(yr,x_m:x_f-1)=tempagb; %[g] adds the year's marsh aboveground biomass to the matrix in its correct row

    x_f = find(elevation(yr-1,:)<(msl(yr)+amp-Dmin),1,'last')+1; %new first forest cell

    %%Update forest soil organic matter
    forestage=forestage+1; %adds a year to forestage to be used in loop 
    df = -msl(yr) + elevation(yr,x_f:B); %depth of forest above mean sea level
    organic_dep_autoch(yr,x_f:B)=f0+fwet*exp(-fgrow*df); %calculate organic deposition in the forest
    mineral_dep(yr,x_f:B)=forestMIN(1,80); %implement mineral deposition in the forest
    aboveground_forest(yr,x_f:B)=BMax_forest./(1+a*exp(-b*df)); %Update forest aboveground biomass
    
    [compaction,tempFd,organic_dep_autoch]=decompose(x_m,x_f,organic_dep_autoch,elevation); %decompose marsh OM, see function decompose for more detail
    Fd(yr)=tempFd; %[kg] Flux of organic matter out of the marsh due to decomposition

    elevation(yr,x_m:B) = elevation(yr,x_m:B)-compaction(x_m:B); %Adjust marsh and forest elevation due to compaction from decomposition
    OM_sum_au(yr,1:length(elevation))=sum(organic_dep_autoch(1:yr,:)); %calculate all auto. OM from the given year
    OM_sum_al(yr,1:length(elevation))=sum(organic_dep_alloch(1:yr,:)); %calculate all allo. OM from the given year

    F=0; %flooding marker
    while x_m < B && x_m<x_f %When the marsh edge is less than the total width of the domain (i.e. the bay-marsh boundary is within the domain)...
        if organic_dep_autoch(yr,x_m) > 0 %If organic deposition is greater than zero, the marsh is still growing
            break
        else %Otherwise, the marsh has drowned, and will be eroded to form new bay
            F=1; %set marker (F) to 1
            edge_flood(yr) = edge_flood(yr) + 1; %Count that cell as a flooded cell
            bfo = bfo + 1;%Increase the bay fetch by one cell
            x_m = x_m + 1;%Update the new location of the marsh edge
        end
    end

    x_f=max(x_m+1,x_f); %forest edge can't be less than or equal to marsh edge 

    if F == 1 %If flooding occurred, adjust marsh flux
        [FF_org,FF_min] = calcFE(bfo,fetch(yr-1),baydepth,elevation,organic_dep_autoch,organic_dep_alloch,mineral_dep,msl,boundyr); %Calculate the amount of organic and mineral sediment liberated from the flooded cells
        Fm_min = Fm_min - FF_min; %Adjust flux of mineral sediment to the marsh
        Fm_org = Fm_org - FF_org; %Adjust flux of organic sediment to the marsh
        elevation(yr,1:x_m) = elevation(yr,1); %Change the drowned marsh cell to z bay cell
    end

    fluxes(:,yr)=[Fe_min Fe_org Fm_min Fm_org Fc_min Fc_org Fb_min Fb_org]'; %[kg/yr] put all fluxes into one matrix for later analysis (saved as model output later)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %Update inputs for marsh edge erosion
    Marsh_edge(yr) = x_m;
    Forest_edge(yr) = x_f;
    Bay_depth(yr) = db;
        
    if x_m > 0 && x_m < B %if the marsh edge is within the domain 
        dmo = msl(yr)+amp-elevation(yr,x_m);%DCW- use following code if you don't want to save elevation through time "elevation(x_m)"
        Edge_ht(yr) = dmo; %save marsh edge elevation in vector
    elseif x_m <= 0 %Condition for if the marsh has expanded to fill the basin.
        disp('Marsh has expanded to fill the basin.')
        endyear = yr;
        edge_flood(endyear+1:end) = []; %empty matrix edge_flood for rest of years
        break
    elseif x_m >= B %Condition for if the marsh has eroded completely away
        disp('Marsh has retreated. Basin is completely flooded.')
        endyear = yr;
        edge_flood(endyear+1:end) = []; %empty matrix edge_flood for rest of years
        break
    end
    
    if db < 0.3 %Condition for if the bay gets very shallow.
        disp('Bay has filled in to form marsh.')
        endyear = yr;
        edge_flood(endyear+1:end) = []; %empty matrix edge_flood for rest of years
        break
    end
    
    Fc_ODE = [];  %empty out the model calculating flux of sediment into and out of bay for next year in for-loop
    C_e_ODE = []; %empty out the edge erosion model for next year in for-loop
    clear tempelevation %clear temporary elevation matrix for next year in for-loop

    fetch(yr)=bfo; %Save change in bay fetch through time
    Fm_flood = 0; %empty out flux from marsh due to flooding

end

marshOM_final=(sum(sum(organic_dep_autoch(:,x_m:x_f-1)))+sum(sum(organic_dep_alloch(:,x_m:x_f-1))))/1000; %[kg] Total mass of organic matter in the marsh at the end of the simulation
marshMM_final=sum(sum(mineral_dep(:,x_m:x_f-1)))/1000; %[kg] Total mass of mineral matter in the marsh at the end of the simulation
marshLOI_final=marshOM_final/(marshOM_final+marshMM_final).*100; %[%] Average loi of the marsh across the marsh platform at the end of the simulation
marshOCP_final=.4.*marshLOI_final+.0025.*marshLOI_final.^2; %Organic carbon content (%) from Craft et al (1991)
marshOC_final=marshOCP_final./100.*(marshOM_final+marshMM_final); %[kg] Mass of organic carbon stored in the marsh at the end of the simulation
save([outputfilename 'marshLOI_final.mat'],'marshLOI_final')
save([outputfilename 'marshOM_final.mat'],'marshOM_final')
save([outputfilename 'marshOC_final.mat'],'marshOC_final')

save([outputfilename 'mortality.mat'],'mortality')
save([outputfilename 'OM_sum_al.mat'],'OM_sum_al')
save([outputfilename 'OM_sum_au.mat'],'OM_sum_au')

forestOM_final=(sum(sum(organic_dep_autoch(:,x_f:B)))+sum(sum(organic_dep_alloch(:,x_f:B))))/1000; %[kg] Total mass of organic matter in the forest at the end of the simulation
forestMM_final=sum(sum(mineral_dep(:,x_f:B)))/1000; %[kg] Total mass of mineral matter in the forest at the end of the simulation
forestLOI_final=forestOM_final./(forestOM_final+forestMM_final).*100; %[%] Average loi of the forest the end of the simulation
forestOCP_final=.4.*forestLOI_final+.0025.*forestLOI_final.^2; %Organic carbon content (%) from Craft et al (1991)
forestOC_final=forestOCP_final./100.*(forestMM_final+forestOM_final); %[kg] total mass of C stored in the forest at the end of the simulation
save([outputfilename 'forestLOI_final.mat'],'forestLOI_final')
save([outputfilename 'forestOM_final.mat'],'forestOM_final')
save([outputfilename 'forestOC_final.mat'],'forestOC_final')

bayexportOM_final=sum(BayExport(:,1)); %[kg] Total mass of organic matter exported from the bay to outside the system over the course of the simulation
bayexportMM_final=sum(BayExport(:,2)); %[kg] Total mass of mineral matter exported from the bay to outside the system over the course of the simulation
bayexportLOI_final=bayexportOM_final./(bayexportOM_final+bayexportMM_final).*100; %[%] Average loi of the material exported from the bay over the course of the simulation
bayexportOCP_final=.4.*bayexportLOI_final+.0025.*bayexportLOI_final.^2; %Organic carbon content (%) from Craft et al (1991)
bayexportOC_final=bayexportOCP_final./100.*(bayexportOM_final+bayexportMM_final); %[kg] total mass of C exported from the bay to outside the system over the course of the simulation
save([outputfilename 'bayexportOM_final.mat'],'bayexportOM_final')
save([outputfilename 'bayexportOC_final.mat'],'bayexportOC_final')
save([outputfilename 'bayexportLOI_final.mat'],'bayexportLOI_final')

bayOM_final=sum(BayOM); %[kg] Total mass of organic matter stored in the bay at the end of the simulation
bayMM_final=sum(BayMM); %[kg] Total mass of mineral matter stored in the bay at the end of the simulation
bayLOI_final=bayOM_final./(bayOM_final+bayMM_final).*100; %[%] Average loi of the material exported from the bay over the course of the simulation
bayOCP_final=.4.*bayLOI_final+.0025.*bayLOI_final.^2; %Organic carbon content (%) from Craft et al (1991)
bayOC_final=bayOCP_final./100.*(bayOM_final+bayMM_final); %[kg] total mass of C stored in the bay at the end of the simulation
save([outputfilename 'bayOM_final.mat'],'bayOM_final')
save([outputfilename 'bayOC_final.mat'],'bayOC_final')
save([outputfilename 'bayLOI_final.mat'],'bayLOI_final')

bgb_final=sum(bgb_sum)/1000; %[kg] Total mass of belowground biomass that has been deposited in the marsh over the course of the simulation
save([outputfilename 'bgb_final.mat'],'bgb_final')

decomp_final = sum(Fd); %[kg] Total mass of organic matter decomposed in the marsh over the course of the simulation
save([outputfilename 'decomp_final.mat'],'decomp_final')

elev=elevation;
save([outputfilename 'OCb.mat'],'OCb')
save([outputfilename 'elev.mat'],'elev')
save([outputfilename 'mineral deposition.mat'],'mineral_dep')
save([outputfilename 'organic deposition.mat'],'organic_dep_autoch','organic_dep_alloch')
save([outputfilename 'marsh edge.mat'],'Marsh_edge')
save([outputfilename 'forest edge.mat'],'Forest_edge')
save([outputfilename 'bay depth.mat'],'Bay_depth')
save([outputfilename 'marsh edge flooding.mat'],'edge_flood')
save([outputfilename 'suspended sediment.mat'],'C_e')
save([outputfilename 'fetch.mat'],'fetch')
save([outputfilename 'compaction.mat'],'compaction')
save([outputfilename 'bgb.mat'],'bgb')
save([outputfilename 'agb.mat'],'agb')
save([outputfilename 'rhomt.mat'],'rhomt')
save([outputfilename 'Fe_org.mat'],'Fe_org')
save([outputfilename 'fluxes.mat'],'fluxes')
save([outputfilename 'aboveground_forest.mat'],'aboveground_forest')
save([outputfilename 'meansealevel.mat'],'msl')
save([outputfilename 'SLR.mat'],'SLR')

%save spinup files that will be imported into transectwithtransectspinup.m function
outputfilefolder=['transectspinupfiles/' filename ];
save([outputfilefolder 'spinupsMarshStrat_initialRSLR' num2str(RSLRi) '_rampedfinalRSLR' num2str(RSLRf) '_CO' num2str(Coi) '.mat'],'msl','SLR','organic_dep_alloch','organic_dep_autoch','mineral_dep','bgb','fluxes','elev') 

clear all