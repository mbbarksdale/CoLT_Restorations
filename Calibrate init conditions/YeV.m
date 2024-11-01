function [Hs,Tp]=YeV(fetch,wind,h)

%Calculates wave height (Hs) and period (Tp) for a given fetch, wind speed, and depth; based on a set of semi-empirical equations 
%from Young and Verhagen (1996)

g=9.8; %Acceleration due to gravity [m/s2]
delta=h*g./wind.^2; %Dimensionless coefficient
chi=fetch*g./wind.^2; %Dimensionless coefficient
epsilon=3.64*10^-3*(tanh(0.493*delta.^0.75).*tanh(3.13*10^-3*chi.^0.57./tanh(0.493*delta.^0.75))).^1.74; %Dimensionless coefficient
ni=0.133*(tanh(0.331*delta.^1.01).*tanh(5.215*10^-4*chi.^0.73./tanh(0.331*delta.^1.01))).^-0.37; %Dimensionless coefficient
Hs=4*sqrt(wind.^4.*epsilon/g^2); %Wave height [m]
Tp=wind./ni/g; %Wave period [s]