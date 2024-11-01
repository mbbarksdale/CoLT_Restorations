function dMWvT(filename)

%plots change in marsh width, as a function of edge erosion and forest retreat, through time for a given RSLR and Co

% close all

% filepath = ['Outputs\' datei '\ws50\RSLR' num2str(RSLR) '/Co' num2str(Co) '_Wind6_150yr'];
outputfilename = ['Run Files\' filename '\Outputs'];
Inputdata=xlsread(['Run Files\' filename '\Input variables.xlsx']);
RSLR=Inputdata(1); %[mm/yr] Relative sea level rise rate
slope=Inputdata(3); %Upland slope
endyear=Inputdata(4); %[yr] Number of years to simulate


% if exist([filepath '/marsh width.mat']) == 0
%     load([filepath '\organic deposition.mat']) %organic_dep - organic deposition in each cell
%     
%     for yr = 2:numel(organic_dep(:,1))
%         MW(yr) = numel(find(organic_dep(yr,:)>0.01));
%     end
%     MW(1) = 1000 + .523/slope;
%     save([filepath '/marsh width.mat'],'marsh width')
% else
    load([outputfilename '/marsh width.mat'])
% end

load([outputfilename '\marsh edge.mat']) %Marsh_edge - organic deposition in each cell
% load([filepath '\marsh edge flooding.mat']) %edge_flood - organic deposition in each cell

FR = (RSLR/1000/slope)*(0:endyear-1);%forest retreat
length(FR)
figure
hold on
minY=min([min(Marsh_edge(1:endyear)-Marsh_edge(1)) min(MW(1:endyear)-MW(1)) min(FR)])-100;
maxY=max([max(Marsh_edge(1:endyear)-Marsh_edge(1)) max(MW(1:endyear)-MW(1)) max(FR)])+100;

set(gca,'Xlim',[0 numel(MW)],'Ylim',[minY maxY])
plot(1:endyear,MW(1:endyear)-MW(1),'-k','LineWidth',2)
plot(1:endyear,Marsh_edge(1:endyear)-Marsh_edge(1),'-g','LineWidth',1)
% plot(1:endyear,cumsum(edge_flood(1:endyear)),'-b','LineWidth',2)
plot(1:endyear,FR,'-r','LineWidth',1)

legend('Marsh Width','Marsh Edge','Forest Edge','Location','NorthWest')

xlabel('Time (yr)','FontSize',15)
ylabel('Change (m)','FontSize',15)

print('-dpng',[outputfilename '\Change in Marsh Width vs Time.png'])