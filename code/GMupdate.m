function [output]=GMupdate(input)
K=input.K;
mu=input.mu;
sig=input.sig;
alp=input.alp;
z=input.z;
N=numel(z);
count=0;
while 1==1
count=count+1;
tem=0;
for j=1:K
        tem=tem+alp(j)*normpdf(z,mu(j),sig(j)); 
end

for j=1:K
        gamz=alp(j)*normpdf(z,mu(j),sig(j))./tem;
        mu(j)=sum(gamz.*(z))/sum(gamz);
        sig(j)=sqrt((gamz.*(z-mu(j)))*(z-mu(j))'/sum(gamz));  
        alp(j)=sum(gamz)/N;
end
clear tem temv
%% delete this part if the target PDF is not symmetric
input.sig=sig;
input.alp=alp;
input.mu=mu;
[output]=EnfSym(input);%% Enforce symmetry of the GM model
sig=output.sig;
alp=output.alp;
mu=output.mu;
%%
pdft=0;
for l=1:K
         pdft=alp(l)*normpdf(z,mu(l),sig(l))+pdft;
end
    ce(count)=-sum(log(pdft))/N;
    if count>=3&&abs((ce(count)-ce(count-1))/ce(count))<1e-8
        break
    end
    if count==K*100
        break
    end
end

output.alp=alp;
output.sig=sig;
output.mu=mu;

end