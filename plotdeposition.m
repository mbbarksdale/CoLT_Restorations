function plotdeposition(filename,ts)

%Function to plot the deposition of mineral and organic matter over the
%marsh surface for a given time step (ts) from a given model run (based on input
%parameters.
    
% filename = ['Outputs/ws' num2str(ws*1000) '/RSLR' num2str(RSLR) '/Co' num2str(Co) '_Wind' num2str(Wind) '_' num2str(Dur) 'yr'];
% filename = ['Outputs/' datein '/ws' num2str(ws*1000) '/RSLR' num2str(RSLR) '/Co' num2str(Co) '_Wind' num2str(Wind) '_' num2str(dur) 'yr'];
outputfilename = ['Run Files\' filename '/Outputs'];
    
% load([filename '/elevation'])
load([outputfilename '/mineral deposition'])
load([outputfilename '/organic deposition'])
load([outputfilename '/accretion'])

L = length(mineral_dep(1,:));
total_acc = total_acc.*1000; %Convert accretion rate from m/yr to mm/yr

% elevation = elevation - (RSLR*ts)/1000; %Set elevation relative to sea level

figure
hold on
subplot(1,3,1)
set(gca,'FontSize',14)
% [hAx,hLine1,hLine2] = plotyy(1:L,elevation(ts,:),1:L,mineral_dep(ts,:));
plot(1:L,mineral_dep(ts,:))
set(gcf,'units','Inches','position',[1 1 18 4.8],'PaperPositionMode','auto')
% set(hAx(1),'ycolor','k')
% set(hAx(2),'ycolor','b')
% set(hLine1,'color','k')
% set(hLine2,'color','b')

%Label the plot/axes
title(['Mineral Deposition at t = ' num2str(ts) ' yr'])
xlabel('Distance (m)')
% ylabel(hAx(1),'Elevation relative to sea level (m)')
% ylabel(hAx(2),'Mineral Deposition (g/m^2)')
ylabel('Mineral Deposition (g/m^2)')

subplot(1,3,2)
set(gca,'FontSize',14)

% [h2Ax,h2Line1,h2Line2] = plotyy(1:L,elevation(ts,:),1:L,organic_dep(ts,:));
plot(1:L,organic_dep(ts,:));
% set(gcf,'units','Inches','position',[7 1 5.8 4.8])
% set(h2Ax(1),'ycolor',[0 0 0])
% set(h2Line1,'color','k')

%Label the plot/axes
title(['Organic Deposition at t = ' num2str(ts) ' yr'])
xlabel('Distance (m)')
% ylabel(h2Ax(1),'Elevation relative to sea level (m)')
% ylabel(h2Ax(2),'Organic Deposition (g/m^2)')
ylabel('Organic Deposition (g/m^2)')

subplot(1,3,3)
set(gca,'FontSize',14)

% [h3Ax,h3Line1,h3Line2] = plotyy(1:L,elevation(ts,:),1:L,total_dep(ts,:));
plot(1:L,total_acc(ts,:));
% set(gcf,'units','Inches','position',[13 1 5.8 4.8])
% set(h3Ax(1),'ycolor','k')
% set(h3Ax(2),'ycolor','r')
% set(h3Line1,'color','k')
% set(h3Line2,'color','r')

%Label the plot/axes
title(['Total Accretion at t = ' num2str(ts) ' yr'])
xlabel('Distance (m)')
% ylabel(h3Ax(1),'Elevation relative to sea level (m)')
% ylabel(h3Ax(2),'Net Accretion (m)')
ylabel('Net Accretion (mm/yr)')
%save the plots
if exist('plottitle') ~= 0
    savefilename = [outputfilename '/' plottitle];
else
    savefilename = [outputfilename '/Deposition plot'];
end

print('-dpng',savefilename)
