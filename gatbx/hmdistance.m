function [m,d] = hmdistance(x)
   n=size(x);
   a=max(x);
for j=1:n 
    b=x(j);
d(j) = sum(xor(a,b));
m(1)=1/n*exp(-d(1));
if j>=2
m(j)=1/n*exp(-d(j))+1/n*exp(-d(1));
end
end

