function [gkk] = autoPSD(w) 
   scal=1/pi;%% scaling factor
%    wfk=15.0;%natural frequency of the soil layer
%    zetafk=0.6;%damping of the soil layer
%    wgk=1.5;%frequency of the second filter
%    zetagk=0.6;%damping of the second filter 
%    gkk=(scal*(wfk^4+4*zetafk^2*wfk^2*w.^2)./((wfk^2-w.^2).^2+4*zetafk^2*wfk^2*w.^2)).*w.^4./((wgk^2-w.^2).^2+4*zetagk^2*wgk^2*w.^2);
gkk=scal*w./w;%% Using the white noise model 
end

