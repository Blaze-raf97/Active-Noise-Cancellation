function [Y,h] = LMSRAF(d,x,N,del)
M = length(x);
y = zeros(1,M);
h = zeros(1,N);
for i = N:1:M
    x1 = x(i:-1:i-N+1);
    y(i) = dot(x1,h);
    e = d(i) - y(i);
    h = h + (del*e).*x1;
end
%h = h./max(h);%normalizing
Y = conv(x,h);
end

