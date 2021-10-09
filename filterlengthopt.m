function [error] = filterlengthopt(x,d)

M = length(x);
NTest = round(M/10):round(2*M/3);%selecting N array for optimization
P = sum(x.^2)/M;
error = zeros(1,length(NTest));
for  i = 1:1:length(NTest)
    del = 1/(10*NTest(i)*P);
    [y,h] = LMSRAF(d,x+xref,NTest(i),del);
    z = y(1:M);
    e = (d-z);
    k = round(M/2)+1;
    error(i) = sum(e(k:M).^2)/(M-k+1);
end

 
end

