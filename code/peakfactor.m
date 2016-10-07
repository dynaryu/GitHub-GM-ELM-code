function [peak] = peakfactor(tk,k) 
nu=tk/2/pi*sqrt(Moment(2,k)/Moment(0,k))/(-log(0.5704));
sigu=sqrt(1-Moment(1,k)^2/Moment(0,k)/Moment(2,k));
peak=sqrt(2*log(2*nu*(1-exp(-sigu^1.2*sqrt(pi*log(2*nu))))));
end