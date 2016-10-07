function [output]=EnfSym(input)
mu=input.mu;
alp=input.alp;
sig=input.sig;
Ko=input.K;

if mod(Ko,2) == 0
[mutem1,ix1]=sort(mu);
ids=[-ones(1,Ko/2),ones(1,Ko/2)];
mutem2=fliplr(mutem1);
mu=(abs(mutem1)+abs(mutem2))/2;
mu=ids.*mu;
alptem1=alp(ix1);
alptem2=fliplr(alptem1);
alp=(alptem1+alptem2)/2/sum((alptem1+alptem2)/2);
sigtem1=sig(ix1);
sigtem2=fliplr(sigtem1);
sig=(sigtem1+sigtem2)/2;

elseif mod(Ko,2)~=0&&Ko~=1
[mutem1,ix1]=sort(mu);
ids=[-ones(1,(Ko-1)/2),0,ones(1,(Ko-1)/2)];
mutem2=fliplr(mutem1);
mu=(abs(mutem1)+abs(mutem2))/2;
mu=ids.*mu;
alptem1=alp(ix1);
alptem2=fliplr(alptem1);
alp=(alptem1+alptem2)/2/sum((alptem1+alptem2)/2);
sigtem1=sig(ix1);
sigtem2=fliplr(sigtem1);
sig=(sigtem1+sigtem2)/2;

else
    mu=0;
end

output.mu=mu;
output.alp=alp;
output.sig=sig;

end