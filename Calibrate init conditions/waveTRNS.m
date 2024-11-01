function W=waveTRNS(amp,depth,D,dist,wind,fetch,hm,hb);

depth=hb;% scarp height
fac=min(1,depth/(2*amp));D=(depth+(depth-fac*2*amp))/2;
[Hs,Tp]=YeV(fetch,wind,D);
kk=wavek(1./Tp,D);cg=2*pi/kk/Tp*0.5*(1+2*kk*D/(sinh(2*kk*D)));
W=cg*9800/16*abs(Hs).^2;