function out = int_exp(x,tau)
%integration of a normalzied exponential, from -inf to x.
out = zeros(size(x));
ind = find(x<0);
if ind
    out(ind) = exp(-abs(x(ind))/tau)/2;
end
ind = find(x>=0);
if ind
    out(ind) = 1 - exp(-abs(x(ind))/tau)/2;
end