function [B,dfo,elevation] = buildtransect(amp,wind,bfo,filename,endyear,R,msl,C,slope,mwo,elev_25,initialspinupRSLR)

spinupdur=width(elev_25'); %determine length of spinup by taking assessing size of the spinup elevation matrix
Dmin=0; %sets minimum depth below MHW at which marsh plants can grow (0 m)
tr=amp*2; %tidal range is 2x amplitude
Dmax=(0.237*tr)-0.092+amp; %maximum depth below MHW at which marsh plants can grow, from McKee and Patrick (1988), flexible for different tidal ranges

fh1 = figure; %create figure to visualize initial domain configuration
hold on

x_m = bfo+1; %marsh edge position is the width of the bay (i.e. the fetch) plus 1
maxY = max(msl) + amp + .5; %gives maximum elevation that the water level will reach by the end of the scenario (assumes a constant rate of SLR) [m]
max_potentialwidth = ceil(maxY/slope); %convert elevation change of water level to the width of forest that will be invaded by water/become marsh [m]
upland_width=ceil(maxY/slope); %automatically calculates the exact amount of forest needed to allow the marsh to migrate 

if max_potentialwidth > upland_width %if you set forest width manually, this checks to make sure the domain is large enough 
    disp('WARNING! Slope/Sea level rise conditions are such that the model domain is too small for the scenario. Adjust the upland width accordingly.')
end

B = bfo + mwo + upland_width; %Total domain width [m], and also number of cells in domain each with 1-m width

x = 1:B; %position of each cell in model domain
elevation = zeros(endyear,B); %create empty matrix to store elevation values
elevation(1:spinupdur,1:x_m+mwo-1) = elev_25(1:spinupdur,1:x_m+mwo-1)-(msl(spinupdur)-msl(1)); %Bay and marsh elevation come from model spin up, adjusted to modern sea level

dfo=msl(spinupdur+1)+amp-elevation(spinupdur,1); %saving bay depth from spinup

%Now form the underlying forest stratigraphy for future organic deposits at depth. 
modernslope = slope*(1:upland_width)+elevation(spinupdur,x_m+mwo-1); %equation for forest slope; takes input forest slope and width of forest, and uses elevation of last marsh cell as the starting elevation of the forest slope

for i = 1:spinupdur %for the length of the spinup
    elevation(i,x_m+mwo:B)=slope*(1:upland_width)+elevation(i,x_m+mwo-1); % implement elevation profile of the forest, always the same slope, using the elevation of the marsh-forest edge to set the starting elevation of the forest
    plot(x,elevation(i,:),'-','Color',[.5 .5 .5]) %plots initial stratigraphy
end

%plots initial stratigraphy and formats figure
plot(x,elevation(spinupdur,:),'k-','LineWidth',1.5) %plots initial stratigraphy
set(gcf,'units','Inches','position',[1 2 12 5],'PaperPositionMode','auto')
set(gca,'FontSize',12)
text(B*.1,max(elevation(spinupdur,:))*.9,['RSLR = ' num2str(R) ' mm/yr'],'FontSize',12)
text(B*.1,max(elevation(spinupdur,:))*.7,['Wind Speed = ' num2str(wind) ' m/s'],'FontSize',12)
text(B*.1,max(elevation(spinupdur,:))*.8,['C = ' num2str(C) ' kg/m^3'],'FontSize',12)
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