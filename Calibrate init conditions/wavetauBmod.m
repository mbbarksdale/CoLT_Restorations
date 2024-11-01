function tw=wavetauBmod(fetch,wind,Dm,Bfrac);
if (Bfrac==0)  ;
[Hs,Tp]=YeV(fetch,wind,Dm);
kk=wavek(1./Tp,Dm);
Um=(pi*Hs./Tp./sinh(kk.*Dm));
aw=Tp*Um/(2*pi);ko=0.001;
fw=0.4*(aw/ko)^-0.75;
tw=1/2*1020*fw*Um^2;
else
tw=0;
end
