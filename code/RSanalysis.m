%% Use peak factor (Vanmarcke's solution) to compute ordinates of the response spectrum 
for k=1:numel(keq)
    peakf=peakfactor(Main.tn,keq(k));
    d(k)=sig(k)*peakf;
    clear y
end
[mzk f1]=sort(d+abs(mu));
alpk=alp(f1);
f2=mzk>=0.95*alp*(d+abs(mu))';
mz=mzk(f2);
maxz=(alpk(f2)/sum(alpk(f2)))*mz'
fprintf('\t Mean peak response:\t%1.4f\n',maxz);