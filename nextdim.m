function dim1 = nextdim(dim)
%find the dimension that a hadamard matrix exist.

d1= min(find(dim./2.^(1:50)<=1));            %the minimum number to satisfy 2.^d > dim.
d2= min(find(dim/12./2.^(0:50)<=1));    %the minimum number to satisfy 12*2.^d > dim.
d3= min(find(dim/20./2.^(0:50)<=1));         %the minimum number to satisfy 20*2.^d > dim.
dim1 = min([2^d1; 12*2^(d2-1); 20*2^(d3-1)]); %increase the dimension as required by hadamard.