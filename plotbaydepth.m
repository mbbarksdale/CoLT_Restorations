function plotbaydepth(filename)

%Plots bay width and depth versus time

% (Bay_width,Bay_depth,endyear,filename)

% filename = ['Outputs/' datein '/ws' num2str(ws*1000) '/RSLR' num2str(RSLR) '/Co' num2str(Co) '_Wind' num2str(Wind) '_' num2str(dur) 'yr'];
outputfilename = ['Run Files\' filename '/Outputs'];

load([outputfilename '/bay depth'])
load([outputfilename '/marsh edge'])
Inputdata=xlsread(['Run Files\' filename '\Input variables.xlsx']);
dur=Inputdata(4); %[yr] Duration of simulation

Bay_width = Marsh_edge - 1;
% Bay_width = Bay_width - 1;
fh2 = figure;
[hAx,hLine1,hLine2] = plotyy(1:dur,Bay_width,1:dur,Bay_depth);
% set(hLine1,'Marker','o')
% set(hLine2,'Marker','x')

%Set the axes limits
space1 = (max(Bay_width) - min(Bay_width))*.25;
uppery1 = max(Bay_width) + space1;
lowery1 = min(Bay_width) - space1/2;
space2 = (max(Bay_depth) - min(Bay_depth))*.25;
uppery2 = max(Bay_depth) + space2;
lowery2 = min(Bay_depth) - space2/2;
ylim(hAx(1),[lowery1 uppery1])
ylim(hAx(2),[lowery2 uppery2])
xlim(hAx(1),[0 dur+.9])
xlim(hAx(2),[0 dur+.9])

%Label the axes
xlabel('Time (yr)')
ylabel(hAx(1),'Bay width (m)')
ylabel(hAx(2),'Bay depth (m)')

%set the tick marks
if uppery1 - lowery1 < 5
    Yincrement1 = 1;
elseif uppery1 - lowery1 < 50
    Yincrement1 = 5;
elseif uppery1 - lowery1 < 100
    Yincrement1 = 10;
elseif uppery1 - lowery1 < 200
    Yincrement1 = 20;
elseif uppery1 - lowery1 < 250
    Yincrement1 = 25;
elseif uppery1 - lowery1 < 500
    Yincrement1 = 50;
elseif uppery1 - lowery1 < 1000
    Yincrement1 = 100;
elseif uppery1 - lowery1 < 2000
    Yincrement1 = 200;
elseif uppery1 - lowery1 < 4000
    Yincrement1 = 400;
else
    Yincrement1 = 1000;
end
% 
Yts1 = Bay_width(1):-Yincrement1:lowery1-1;
Yticks1 = Yts1(end):Yincrement1:uppery1+1;
set(hAx(1),'YTick',Yticks1)
% set(hAx(1),'YLim',[1300 5200],'YTick',1500:500:5000)%Set to match M & C expansion plots
% set(hAx(1),'XLim',[0 5000],'YLim',[4800 10200],'YTick',5000:500:10000)%Set to match M & C retreat plots

if uppery2 - lowery2 < .1
    Yincrement2 = .01;
elseif uppery2 - lowery2 < .2
    Yincrement2 = .02;
elseif uppery2 - lowery2 < .5
    Yincrement2 = .05;
elseif uppery2 - lowery2 < 1
    Yincrement2 = .1;
elseif uppery2 - lowery2 < 2
    Yincrement2 = .2;
elseif uppery2 - lowery2 < 3
    Yincrement2 = .25;
else
    Yincrement2 = .5;
end

Yts2 = Bay_depth(1):-Yincrement2:lowery2-1;
Yticks2 = Yts2(end):Yincrement2:uppery2+1;
set(hAx(2),'YTick',Yticks2)
% set(hAx(2),'YLim',[1.45 2.25],'YTick',1.5:.1:2.2) %Set to match M & C expansion plots
% set(hAx(2),'XLim',[0 5000],'YLim',[1.75 3.05],'YTick',1.8:.1:3)%Set to match M & C retreat plots

%print the final width and depths on the plot
text(.1*dur,max(Bay_width)+space1*2/3,['Final bay width = ' num2str(Bay_width(end)) ' m'])
text(.1*dur,max(Bay_width)+space1*1/3,['Final bay depth = ' num2str(Bay_depth(end)) ' m'])
legend('Bay Width','Bay Depth')
set(gcf,'units','Inches','position',[1 1 5.8 4.8])

%save the plots
savefilename = [outputfilename '/Bay depth & width plot'];
print('-dpng',savefilename)
saveas(fh2,[savefilename '.fig'])