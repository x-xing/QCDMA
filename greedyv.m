function out = greedyv(v,dim)
%greedy algorithm to find the set of v that is maximum for xor.

%generate the goal function. 
ip = floor((length(v)-1)/2)+1; 
out = ip*(ip-rem(length(v),2))*(dim-1);
