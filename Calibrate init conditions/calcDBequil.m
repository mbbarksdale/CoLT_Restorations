function db_eq = calcDBequil

%Function to find the initial condition for bay depth at which the model
%will experience the lowest rate of change. Outputs an array of bay depths
%for a given combination of rate of sea level rise and external sediment
%supply. Will produce a different array for each fetch and wind speed
%conditions.

close all

fetch = 5000; %mudeflat width [m]
wind = 6; %reference wind speed [m/s]

filename = ['Fetch' num2str(fetch) '_Wind' num2str(wind)];

if exist(filename) == 0
    mkdir(filename)
end

fh1 = figure;
set(fh1,'units','inches','position',[0.1,.45,19.8,9.9])

for R = 1:15
    subplot(3,5,R)
    hold on
    title(['RSLR = ' num2str(R) 'mm/yr'])
    for C = 10:10:150
        [R C]
        [DM DB] = wetland3P(R,C,wind,fetch);
        plot(1:numel(DB),DB,'k-')
        xlabel('Time (yr)')
        ylabel('Bay Depth')
        for i = 1:numel(DB)-10
                if sign(DB(i+9)-DB(i+8)) == sign(DB(i+1)-DB(i))
                    p = polyfit((i:i+9),DB(i:i+9)',1);
                    m(i) = p(1);
                    b(i) = p(2);
                else
                    m(i) = 100;
                    b(i) = 100;
                end
        end
        dbi = find(abs(m)==min(abs(m)));
        db(R,C/10) = b(dbi);
        plot([0 200],[b(dbi) 200*m(dbi)+b(dbi)],'bx-')
        clear dbi m b p        
    end
end

saveas(fh1,[filename '/Change in Bay Depth.fig'])
print('-dpng',fh1,[filename '/Change in Bay Depth.png'])

fh2 = figure;
set(fh2,'units','inches','position',[0.1,.45,19.8,9.9])

subplot(3,5,1)
hold on
for C = 1:15
    p= polyfit(1:10,db(1:10,C)',1);
    db_eq(C,1:15)=polyval(p,1:15);
    subplot(3,5,C)
    plot(1:15,db(:,C))
    hold on
    plot(1:15,db_eq(C,1:15),'r-')
    xlabel('RSLR (mm/yr)')
    ylabel('Equil. Depth (m)')
    ylim([1.75 2.65])
    text(5,1.86,['C_o = ' num2str(C*10)])
end

saveas(fh2,[filename '/db regression.fig'])
print('-dpng',fh2,[filename '/db regression.png'])

fh3 = figure;
xlabel('RSLR (mm/yr)','FontSize',15)
ylabel('Co (kg/m^3)','FontSize',15)
hold on
surf(db_eq)
set(gca,'YLim',[1 15],'XLim',[1 15],'YTick',1:15,'YTickLabel',10:10:150,'XTick',1:15)
colorbar

saveas(fh3,[filename '/DB_phase.fig'])
print('-dpng',fh3,[filename '/DB_phase.png'])

save([filename '/Equilibrium Bay Depth.mat'],'db_eq')
