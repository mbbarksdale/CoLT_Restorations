function [compaction,Fd,organic_dep_autoch] = decompose(x_m,x_f,organic_dep_autoch,elevation)

%Decomposes all of the organic sediment within the marsh soil profile at a rate determined by depth.

global yr
global B
global mui
global mki
global rhoo

compaction = zeros(1,B); %create empty vector for compaction
Fd=0; %set decomp flux to zero

%Decompose the marsh sediment
for x = x_m:x_f-1 %Loop through each marsh and upland cell in the domain
    for tempyr = yr:-1:1 %Loop through each pocket of sediment in each cell, starting at the most recently deposited packet of sediment at the surface
        depth = elevation(yr,x) - elevation(tempyr,x); %Depth of sediment pocket below the surface
        if depth > mui %If depth of sediment pocket is greater than the maximum depth at which decomposition occurs...
            decomp(tempyr)=0; %...then decomp goes to 0...
            break %...and loop breaks
        else %If depth of sed pocket is equal to or less than max depth at which decomp occurs...
            decomp(tempyr) = organic_dep_autoch(tempyr,x)*(mki * exp(-depth/mui)); %[g] Mass of organic material decomposed  from a given "pocket" of sediment
            organic_dep_autoch(tempyr,x) = organic_dep_autoch(tempyr,x) - decomp(tempyr); %[g] Autochthonous organic material in a given "pocket" of sediment updated for deomposition
        end
    end
    compaction(x)=sum(decomp)/1000/rhoo; %[m] Total compaction in a given cell is a result of the sum of all decomposition in that cell
    Fd=Fd+sum(decomp); %[kg] Flux of organic matter out of the marsh due to decomposition
    clear decomp
end