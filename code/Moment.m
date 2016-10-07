function [output] = Moment(n,k)
global m c sc dw Nf
w=dw:dw:dw*Nf;
w=w-dw/2;
output=0;
for p=1:Nf
 H(p)=sc/(k+c*w(p)*1i-m*w(p)^2);
 if p==1
 output=output+w(p)^n*autoPSD(w(p))*dw*abs(H(p))^2;
 else
 output=output+0.5*(w(p)^n*autoPSD(w(p))*dw*abs(H(p))^2+w(p-1)^n*autoPSD(w(p-1))*dw*abs(H(p-1))^2);  
 end
end
output=output*2;
end