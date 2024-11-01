function plotsurface(filename,ts)

% close all

%Function to plot the surface profile for a given time step (ts) from a given model run (based on input
%parameters.

incr = 10; %Increment between plotted surface profiles

if nargin == 0
    RSLR = input('RSLR = ');
    Co = input('Co = ');
    Wind = input('Wind = ');
    Dur = input('Duration = ');
    ts = input('Time Step = ');
    plottitle = input('Plot Title = ');
end

% filename = ['Outputs/RSLR' num2str(RSLR) '/Co' num2str(Co) '_Wind' num2str(Wind) '_' num2str(Dur) 'yr'];
% filename = ['Outputs/' direct '/ws50/RSLR' num2str(RSLR) '/Co' num2str(Co) '_Wind6_150yr'];
% filename = ['Outputs/Paper/Figure 3d sensitivity analysis/WM_' num2str(Wm) '__Wf_' num2str(Wf) '/RSLR_' num2str(RSLR)]; %For sensitivity analysis runs
outputfilename = ['Run Files\' filename '/Outputs'];

load([outputfilename '/elevation'])
load([outputfilename '/marsh edge'])
%Load in input variables
Inputdata=xlsread(['Run Files\' filename '\Input variables.xlsx']);
RSLRi=Inputdata(1);
Coi=Inputdata(2);
slope=Inputdata(3); %Upland slope
endyear=Inputdata(4); %[yr] Number of years to simulate

L = length(elevation(1,:));
x = 1:L;

x = x-5000;
% elevation = elevation - (RSLR*ts)/1000; %Set elevation relative to sea level

figure
hold on

% upperY = max(max(elevation(1:ts,:)));
% lowerY = min(min(elevation(1:ts,:)));
% Ydiff = upperY - lowerY;

% text(L*.4,upperY-Ydiff*.05,['RSLR = ' num2str(RSLR) ' mm/yr'],'FontSize',12)
% text(L*.1,upperY-Ydiff*.1,['Wind Speed = ' num2str(Wind) ' m/s'],'FontSize',12)
% text(L*.4,upperY-Ydiff*.16,['C_o = ' num2str(Co) ' mg/L'],'FontSize',12)

plottitle = ['Time Step ' num2str(ts)];
title(plottitle)

% text(2200,.5,['SSC = ' num2str(Co) ' mg/L'],'FontSize',14,'FontName','Calibri (body)')
% text(2200,0,['SLR = ' num2str(RSLR) ' mm/yr'],'FontSize',14,'FontName','Calibri (body)')
% text(0,2,[plottitle],'FontSize',16,'FontName','Calibri (body)','FontWeight','bold')

set(gcf,'units','Inches','position',[1 1 3.5*2 2*2],'PaperPositionMode','auto')

% plot(x,elevation(1,:),'--k')
plot(x(Marsh_edge(1)),elevation(1,Marsh_edge(1)),'xg','MarkerSize',8)
plot(x(find(elevation(1,:)>=.7,1,'first')),elevation(1,find(elevation(1,:)>=.7,1,'first')),'xr','MarkerSize',8)

h = plot([-1 -1],[-5 -5],'k--');

hl = legend('marsh edge','upland edge',['initial surface \color{white}l'],'Location','SouthEast');
hlc = get(hl,'children');
XD = get(hlc(2),'XData');
set(hlc(2),'XData',[XD(1) XD(1)+(XD(2)-XD(1))*.76])

h2 = plot([-5000 0],[0 0],'b--');
h3 = plot([-5000 x(Marsh_edge(ts))],[RSLRi*ts/1000 RSLRi*ts/1000],'b-');
% h3 = plot([-5000 3500],[RSLR*ts/1000 RSLR*ts/1000],'b-');

% delete(h)

for yr = 1+incr:incr:ts
    plot(x,elevation(yr,:),'Color',[1-yr/(ts+incr) 1-yr/(ts+incr) 1-yr/(ts+incr)]) %For color gradient
%     plot(x,elevation(yr,:),'Color','k')
end

for yr = 1+incr:incr:ts
    plot(x(Marsh_edge(yr)),elevation(yr,Marsh_edge(yr)),'xg','MarkerSize',8)
    plot(x(find(elevation(yr,:)>=RSLRi/1000*yr+.7,1,'first')),elevation(yr,find(elevation(yr,:)>=RSLRi/1000*yr+.7,1,'first')),'xr','MarkerSize',8)
end

plot(x,elevation(ts,:),'-k')
if Marsh_edge(ts) <= L
    plot(x(Marsh_edge(ts)),elevation(ts,Marsh_edge(ts)),'xg','MarkerSize',8);
    plot(x(find(elevation(ts,:)>=RSLRi/1000*ts+.7,1,'first')),elevation(ts,find(elevation(ts,:)>=RSLRi/1000*ts+.7,1,'first')),'xr','MarkerSize',8)
end

plot(x,elevation(1,:),'--k')

colormap gray
% h = colorbar('YTick',[65 1],'YTickLabel',{'Year 1',['Year ' num2str(ts)]});
set(gca,'FontSize',14,'FontName','Calibri (body)','XTick',-4000:1000:3000,'YTick',-4:4)

%Label the plot/axes
% title(['Profile Surface at t = ' num2str(ts) ' yr'])
xlabel('Distance from initial marsh edge (m)','FontSize',16,'FontName','Calibri (body)','FontWeight','bold')
ylabel('Elevation (m)','FontSize',16,'FontName','Calibri (body)','FontWeight','bold')
xlim([-500 3750])
ylim([-2 2.5])

h3 = plot([-5000 x(Marsh_edge(ts))],[RSLRi*ts/1000 RSLRi*ts/1000],'b-');
% h3 = plot([-5000 3500],[RSLR*ts/1000 RSLR*ts/1000],'b-');

print('-dpng',[outputfilename '/Surface Plot ' num2str(ts)])

% print('-dpng',['Plots/Figure 2/' plottitle])
% print('-dmeta',['Plots/' plottitle])
% saveas(gcf,['Plots/AGU Poster/Surface plots/' plottitle '.fig'])

% set(gcf,'units','Inches','position',[1 1 9 7.5],'PaperPositionMode','auto')
% ch = colorbar;
% set(ch,'YTick',1:64/6:65,'YTickLabel',[150 125 100 75 50 25 1])
% print('-dpng',['Plots/Colorbar'])