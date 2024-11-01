function cfphase(filename,a,b,c)

for i=a:b:c
outputfilename = ['Run Files\' filename '\Outputs_RSLR7_CO' num2str(i) '_Slope005'];
Inputdata=xlsread(['Run Files\' filename '\Input variables.xlsx']);

load ([outputfilename '\bayexportOC_final.mat'])

A(7,i/10)=bayexportOC_final;
end

for i=a:b:c
outputfilename = ['Run Files\' filename '\Outputs_RSLR3_CO' num2str(i) '_Slope005'];
Inputdata=xlsread(['Run Files\' filename '\Input variables.xlsx']);

load ([outputfilename '\bayexportOC_final.mat'])

A(3,i/10)=bayexportOC_final;
end

for i=a:b:c
outputfilename = ['Run Files\' filename '\Outputs_RSLR10_CO' num2str(i) '_Slope005'];
Inputdata=xlsread(['Run Files\' filename '\Input variables.xlsx']);

load ([outputfilename '\bayexportOC_final.mat'])

A(10,i/10)=bayexportOC_final;
end

A

figure
pcolor(A)

figure
plot(A(3,:))
hold on
plot(A(7,:))
hold on
plot(A(10,:))

end

