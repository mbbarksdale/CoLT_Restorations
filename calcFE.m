function [FE_org,FE_min] = calcFE(bfoc,bfop,baydepth,elevation,organic_dep_autoch,organic_dep_alloch,mineral_dep,msl,boundyr)

global yr
global rhos
global amp

%Function to calculate the flux of organic matter (FE_org) and the flux of mineral sediment (FE_min) from the marsh to the bay,using
%the fetch for the current year (bfoc) the fetch for the previous year (bfop) and the stratigraphy of organic and mineral deposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Calculate OM eroded from the marsh platform
organic_dep = organic_dep_autoch + organic_dep_alloch; %calculates OM deposition across the marsh as sum of alloch and autoch OM
pyr = yr-1; % defines the previous year
E=bfoc-bfop; %Amount of erosion between the previous yr and the current yr

x_m1 = ceil(bfop); %first marsh cell to erode
x_m2 = ceil(bfoc); %last marsh cell to erode

if  E <= 0 %If there is no erosion...
    FE_org = 0; %...then no OM is eroded [g] 
    FE_min = 0; %...then no mineral matter is eroded [g]
elseif x_m1 == x_m2 %If actively eroding marsh edge has not changed since previous year
    if  boundyr==1 %If the first-year marsh sediment cohort under the current marsh edge is positioned above the previous year's bay bottom --> then would erode entirely through the marsh and into underlying strat 
        us=elevation(1,x_m1)-elevation(yr-1,1); %[m] Depth of underlying stratigraphy
        usmass=us*rhos*1000; %[g] Mass of pure mineral sediment underlying marsh at marsh edge
        FE_org = sum(organic_dep(1:pyr,x_m1))*E; %[g] OM eroded is equal to the total amount of OM in the eroding marsh edge (including both initial deposit and OM deposited since the model run began) times the fraction of the marsh edge cell that is eroded
        FE_min = sum(mineral_dep(1:pyr,x_m1))*E+usmass; %[g] MIN eroded is equal to the total amount of OM in the eroding marsh edge (including both initial deposit and OM deposited since the model run began) times the fraction of the marsh edge cell that is eroded
    else %If depth of erosion is less than marsh deposits
        FE_org = sum(organic_dep(boundyr:pyr,x_m1))*E; %[g] Total mass of OM deposited in the marsh edge cell
        FE_min = sum(mineral_dep(boundyr:pyr,x_m1))*E; %[g] Total mass of MIN deposited in the marsh edge cell
    end
else %If actively eroding marsh edge has changed since previous year
    Hfrac_ero=ones(1,x_m2); %create vector
    Hfrac_ero(1) = x_m1 - bfop; % Horizontal fraction of previous marsh edge that is eroded
    Hfrac_ero(end) = bfoc-floor(bfoc); %Horizontal fraction of current marsh edge that is eroded    
    FE_org=0; %reset FE_org to 0 b/c updated below
    FE_min=0; %reset FE_min to 0 b/c updated below
    for x_m = x_m1:x_m2 %start a for-loop starting with first marsh cell to last marsh cell to erode 
        if  boundyr==1 %If the first-year marsh sediment cohort under the current marsh edge is positioned above the previous year's bay bottom --> then would erode entirely through the marsh and into underlying strat 
            us=elevation(1,x_m)-elevation(yr-1,1); %[m] Depth of underlying stratigraphy
            usmass=us*rhos*1000; %[g] Mass of pure mineral sediment underlying marsh at marsh edge
            FE_org = FE_org + sum(organic_dep(1:pyr,x_m))*Hfrac_ero(x_m); %[g] OM eroded from previous marsh edge cell
            FE_min = FE_min + sum(mineral_dep(1:pyr,x_m))*Hfrac_ero(x_m) + usmass; %[g] MIN eroded from previous marsh edge cell
        else %If depth of erosion is less than marsh deposit (don't have to account for underlying strat in erosion mass) 
            FE_org = FE_org + sum(organic_dep(boundyr:pyr,x_m))*Hfrac_ero(x_m); %[g] OM eroded from previous marsh edge cell
            FE_min = FE_min + sum(mineral_dep(boundyr:pyr,x_m))*Hfrac_ero(x_m); %[g] MIN eroded from the marsh edge cell
        end    
    end
end