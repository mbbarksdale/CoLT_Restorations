function [B,dfo,elevation] = buildtransectspinup(amp,wind,bfo,filename,endyear,R,msl,C,slope,mwo,tempelevation)

%function to establish the intial transect elevation of year 1 of the spinup

spinupdur=1; %misnomer; basically sets the first year of the spinup as the imported data for building the transect
Dmin=0; %defines minimum depth below MHW at which marsh plants can grow (0 m)
tr=amp*2; %defines tidal range based on amp which is imported
Dmax=(0.237*tr)-0.092+amp; %maximum depth below MHW at which marsh plants can grow, from McKee and Patrick (1988), flexible for different tidal ranges

fh1 = figure; %create figure to visualize initial domain configuration
hold on

filename1 = ['Calibrate init conditions/Fetch' num2str(bfo) '_Wind' num2str(wind)]; %name of file with initial conditions to load (in subsequent lines)

%Determine initial bay depth, such that change in depth will be small
if exist(filename1) == 0 %if initial conditions file doesn't exist...
    warning('Initial conditions have not been calibrated for these parameters') %display this error message
    dfo = 2; %set bay depth to 2 m
else %if there is an initial conditions file, load those
    load([filename1 '/Equilibrium Bay Depth.mat']) %load bay depth vector
    load([filename1 '/Equilibrium Marsh Depth.mat']) %load marsh depth vector
    if C/10 > size(db_eq,1) || R > size(db_eq,2)  %if requested model run is outside the range of initial conditions that were calibrated (10-100 mg/L or SLR 1-15 mm/yr)
        warning('Initial conditions have not been calibrated for these parameters 1') %display this error message
        dfo = 2; %set bay depth equal to 2 
    elseif C/10 < .5 || R < .5 %if SSC is less than 5 mg/L or SLR is less than 0.5 mm/yr
        warning('Initial conditions have not been calibrated for these parameters 2') %display this error message
        dfo = 2; %set bay depth equal to 2 m
    else
        dfo = db_eq(round(C/10),round(R)); %rounds both SSC and SLR inputs to integer, then uses lookup table to set the bay depth
    end
end

x_m = bfo+1; %marsh edge position is the width of the bay (i.e. the fetch) plus 1

maxY = max(msl) + amp + .5; %gives maximum elevation that the water level will reach by the end of the scenario (assumes a constant rate of SLR) [m]
max_potentialwidth = ceil(maxY/slope); %convert elevation change of water level to the width of forest that will be invaded by water/become marsh [m]
upland_width=ceil(maxY/slope); %k added, this is what it was in the previous version

if max_potentialwidth > upland_width %if you set forest width manually, this checks to make sure the domain is large enough
    disp('WARNING! Slope/Sea level rise conditions are such that the model domain is too small for the scenario. Adjust the upland width accordingly.') %display this error message
end

B = bfo + mwo + upland_width; %Total domain width [m], and also number of cells in domain each with 1-m width

x = 1:B; %x=position of each cell in model domain
elevation = zeros(endyear,B); %create empty matrix to store elevation values
elevation(1,1:x_m-1) = amp-dfo; %Bay depth for spinup duration is at equilibrium
m=.001; %initial marsh slope
elevation(1,x_m:x_m+mwo-1) = (1:500).*m+amp+.02-Dmax; %builds initial stratigraphy of 500-m wide marsh by setting initial slope and y-intercept

%Now form the underlying forest stratigraphy for future organic deposits at
modernslope = slope*(1:upland_width)+elevation(spinupdur,x_m+mwo-1); %equation for forest slope; takes input forest slope for width of forest, and uses elevation of last marsh cell as y-intercpt
elevation(1,x_m+mwo:B) = modernslope; %-((spinupdur-i)*.025);%(spinupdur/1000));

%create figure to visualize domain configuration
plot(x,elevation(spinupdur,:),'k-','LineWidth',1.5) %plots initial stratigraphy
set(gcf,'units','Inches','position',[1 2 12 5],'PaperPositionMode','auto')
set(gca,'FontSize',12)
text(B*.1,max(elevation(1,:))*.3,['d_b = ' num2str(round(100*dfo)/100) ' m'],'FontSize',12)
text(bfo*.4,0,'Fetch')
text(bfo*.4,-.2,[num2str(bfo) 'm'])
text(bfo+mwo*.2,amp+.25,'Platform')
text(bfo+mwo*.2,amp+.05,[num2str(mwo) 'm'])
text(bfo+mwo+upland_width*.3,elevation(spinupdur,round(bfo+mwo+upland_width*.3))+.8,'Upland Slope')
text(bfo+mwo+upland_width*.3,elevation(spinupdur,round(bfo+mwo+upland_width*.3))+.6,[num2str(slope)])
plot([x(x_m) x(x_m+mwo) x(x_m+mwo)],[elevation(spinupdur,x_m) elevation(spinupdur,x_m+mwo) elevation(spinupdur,x_m+mwo)],'xg','MarkerSize',10,'LineWidth',2)
ylabel('Elevation Relative to Initial Sea Level (m)')
xlabel('Distance (m)')

outputfilename = [filename 'Initial Surface']; %saves output graph in designated folder