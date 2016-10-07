function [sam] = gE(n)
samp=[];
while 1==1
x=-3+6*rand(n,1);
y=rand(n,1);
id=(y<=exp(-x.^4/4)/2.5637);
samp=[samp;x(id)];
if size(samp,1)>=n
    break
end    
end
sam=samp(1:n);
end