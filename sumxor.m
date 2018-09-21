function out = sumxor(v)
%calculates the sum of XORs for nchoosek(v,2).
c=nchoosek(v,2);
out = sum(bitxor(c(:,1),c(:,2)));