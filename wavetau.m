function [tw,Um,Hs,kk] = wavetau(fetch,wind,Df)

%Calculates wave bed shear stress, as a function of fetch, wind speed, and
%bay depth. From Mariotti and Fagherazzi (2013)

[Hs,Tp] = YeV(fetch,wind,Df); %Significant wave height (Hs, [m]) and wave period (Tp, [m]
kk = wavek(1./Tp,Df);%Calculates wave number [m^-1]
Um = (pi*Hs./Tp./sinh(kk.*Df)); %Term in equation for shear stress [m/s]
aw = Tp*Um/(2*pi);%Term in equation for shear stress [m]
ko=0.001; %Roughness [m] (Mariotti and Fagherazzi, 2013)
fw = 0.4*(aw/ko)^-0.75; %Friction factor, dimensionless
tw = 1/2*1020*fw*Um^2; %Wave bed shear stress [Pa]