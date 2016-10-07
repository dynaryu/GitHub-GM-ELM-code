global sigk m c sc dw Nf
dw=Main.dw;
Nf=Main.Nf;
sc=1;%% see Eq.(18) in the GM-ELM paper
m=1; %% see Eq.(18) in the GM-ELM paper
c=1; %% see Eq.(18) in the GM-ELM paper
dth=100;
thr=3/dth:3/dth:3;%%define the thresholds at which the first-passage probabilities are of interest
fp=zeros(1,numel(thr));
x0=[1];
for l=1:Ko
sigk=sig(l);
keq(l)=fminsearch(@objectf,x0);
end
keq=abs(keq);
sigk=siggm;
kgm=fminsearch(@objectf,x0);

for tr=1:size(thr,2)
    z0=thr(tr);
    tem=0;
for kk=1:Ko
    lam0=Moment(0,keq(kk));%% Compute the spectral moment
    lam1=Moment(1,keq(kk));
    lam2=Moment(2,keq(kk));
    tem(kk)=1/(2*pi)*sqrt(lam2/lam0)*exp(-0.5*(z0-mu(kk))^2/lam0); %% Compute the crossing rate for each sublinear system
end
cr(tr)=alp*tem'; %% compute the crossing rate  
end 
fpc=1-exp(-cr*Main.tn);%% compute the first passage probability
thr=thr/sqrt(0.676);
tcr = load(['CR_MCS.txt']);%%load the MCS solution for crossing rate
ecr=load(['CR_ELM.txt']);%%load the ELM solution for crossing rate
figure (4)
semilogy(thr,(ecr),'-sr',thr,(cr),'-o',thr,(tcr),'k','LineWidth',1)
ylim([1e-4 1])
xlim([0 3.6])
legend({'ELM','GM-ELM','MCS'},'FontSize',12)
xlabel('z/\sigma','FontSize',12)
ylabel('Crossing rate','FontSize',12)
tfpc = load(['FP_MCS.txt']);%%load the MCS solution for FP probability
efpc=load(['FP_ELM.txt']);%%load the ELM solution for FP probability
figure (5)
semilogy(thr,(efpc),'-sr',thr,(fpc),'-o',thr,(tfpc),'k','LineWidth',1)
ylim([1e-4 1])
xlim([0 3.6])
legend({'ELM','GM-ELM','MCS'},'FontSize',12)
xlabel('z/\sigma','FontSize',12)
ylabel('First-passage probability','FontSize',12)