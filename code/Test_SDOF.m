%% This code uses GM-ELM to solve the 1st example in Wang Z.& Song J. 2016 Structural Safety
%% For other problems, one needs to modify the 'gE.m' function for response samples generating,
%% and other relevant parameter settings
clear all
clc
Main.tn=35; %% duration of the input stochastic process
Main.cut=5*pi;%% cut-off frequency of the process
Main.dw=Main.cut/200;%% frequency step size
Main.Nf=round(Main.cut/Main.dw);%% number of discretized frequency points
AnalyzeMode=1;%%=1, set K=20; else searching for optimal K starting from K=1
zf=gE(1.0e5)';%% generate response sample points via the analytical PDF solution
bd=std(zf);
%% Start the maximum likelihood estimation process
if AnalyzeMode==1
K=20;
input.alp=1/K*ones(1,K);
input.sig=bd*ones(1,K);
input.mu=-bd:2*bd/(K-1):bd;
input.z=zf;
input.K=K;
output=GMupdate(input);
input.sig=output.sig;
input.alp=output.alp;
input.mu=output.mu;
MU{K}=input.mu;
ALP{K}=input.alp;
SIG{K}=input.sig;
Ko=K;
else
K=0;
count=0;     
while K<=50
    K=K+1
    count=count+1;
    input.K=K;    
    if K==1
     input.alp=1/K*ones(1,K);
     input.sig=bd*ones(1,K);
     input.mu=0; 
    else
     input.alp=1/K*ones(1,K);
     input.sig=bd*ones(1,K);
     input.mu=-bd:2*bd/(K-1):bd; 
    end
input.z=zf;
output=GMupdate(input);
input.sig=output.sig;
input.alp=output.alp;
input.mu=output.mu;
     pdft=0;
     for l=1:K
         pdft=input.alp(l)*normpdf(zf,input.mu(l),input.sig(l))+pdft;
     end
    ce(K)=-sum(log(pdft))/numel(zf);
    MU{K}=input.mu;
    ALP{K}=input.alp;
    SIG{K}=input.sig;
    if count>=3
    if ce(K-1)<ce(K-2)&&ce(K-1)<ce(K)
        Ko=K-1;
        break
    elseif ce(K-2)<ce(K-1)&&ce(K-2)<ce(K)
        Ko=K-2;
        break
    elseif abs(ce(K)-ce(K-1))/abs(ce(K))<1e-6
        Ko=K;
        break
    end
    end
end
end

mu=MU{Ko};
alp=ALP{Ko};
sig=SIG{Ko};
mmu=alp*mu';

fprintf('\t optimal K:\t%1.0f\n', Ko);

for j=1:300
    pdf(j)=sum(alp.*normpdf(-3+0.02*j,mu,sig));
end
xx=-3+.02:0.02:3;
tru=@(x)exp(-x.^4/4);
for j=1:numel(xx)
    cdf(j)=sum(alp.*normcdf(xx(j),mu,sig));
    cdft(j)=integral(tru,-inf,xx(j))/integral(tru,-inf,inf);
end
figure (2)
plot(xx,pdf,'-.b',xx,normpdf(xx,0,sqrt(0.5774)),'--g',xx,exp(-xx.^4/4)/integral(tru,-inf,inf),'r',xx,normpdf(xx,0,(0.8222)),':k','LineWidth',2)
legend({'GM-ELM','ELM','Exact','Gaussian'},'FontSize',12)
xlabel('Displacement','FontSize',12)
ylabel('Probability density','FontSize',12)
figure (3)
semilogy(xx,(1-cdf),'-.b',xx,(1-normcdf(xx,0,sqrt(0.5774))),'--g',xx,(1-cdft),'r',xx,(1-normcdf(xx,0,(0.8222))),':k','LineWidth',2)
xlabel('Displacement','FontSize',12)
ylabel('Complementary cumulative distribution','FontSize',12)
legend({'GM-ELM','ELM','Exact','Gaussian'},'FontSize',12)
vargm=sum(alp.*sig.^2+alp.*(mu-mmu).^2);
siggm=sqrt(vargm);
fprintf('\t Variance:\t%1.4f\n', vargm);
fprintf('\t Standard deviation:\t%1.4f\n', siggm);
run RVanalysis %%run random vibration analysis
run RSanalysis %%run response spectrum analysis
% figure (5)
% plot(-maxz:maxz/1000:maxz,(-maxz:maxz/1000:maxz).^3,'LineWidth',2)
% xlabel('Deformation','FontSize',12)
% ylabel('Restoring force','FontSize',12)