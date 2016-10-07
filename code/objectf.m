function [out] = objectf(k)
global sigk m c sc dw Nf
w=dw:dw:dw*Nf;
w=w-dw/2;
tem=0;

for p=1:numel(w)
H(p)=sc/(k-m*w(p)^2+c*1i*w(p));
 if p==1
 tem=tem+autoPSD(w(p))*dw*abs(H(p))^2;
 else
 tem=tem+0.5*dw*(autoPSD(w(p))*abs(H(p))^2+autoPSD(w(p-1))*abs(H(p-1))^2);  
 end  
end
tem=tem*2;
out=abs(tem-sigk^2)/sigk^2;
