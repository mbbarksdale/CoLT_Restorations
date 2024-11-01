function W = waveTRNS(amp,wind,fetch,hb)

%Calculates wave power density at marsh boundary. From Mariotti and Carr (2014)

depth = hb;% scarp height
fac = min(1,depth/(2*amp));%Proportion of tide that the bay is flooded, dimensionless
D=(depth+(depth-fac*2*amp))/2;%[m] average bay depth over tidal cycle
[Hs,Tp] = YeV(fetch,wind,D); %Solves for wave height (Hs, [m]) and wave period (Tp, [s])
kk = wavek(1./Tp,D); %Solves for wave number [m^-1]
cg=2*pi/kk/Tp*0.5*(1+2*kk*D/(sinh(2*kk*D)));%Wave group celerity at marsh edge [m/s]
W = cg*9800/16*abs(Hs).^2;%[W] Wave power density at the marsh boundary