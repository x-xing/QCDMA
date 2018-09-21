%%XOR of exhaustive search
disp('optimized:');
vxor = zeros(15,1);
com1 = zeros(15,1);
for i1=1:15
    ind1 = indcell{i1}-1;
    linestr = [];
    for i2=1:i1+1
        linestr = [linestr dec2bin(ind1(i2),4) ' '];
    end
    %disp(linestr);
    disp(indcell{i1}-1);
    vxor(i1) = sumxor(ind1);
    com1(i1) = nchoosek(i1+1,2);
end
disp('Sum of XOR:');disp(vxor);
disp('averaged XOR:');disp(vxor./com1);  %result of bitwise XOR.

%% Approximation
% c=nchoosek(0:15,2);
% n1=zeros(1,length(c));
% for i1=1:length(c)
%     c1=nchoosek(c(i1,:),2);
%     n1(i1)= sum(bitxor(c1(:,1),c1(:,2))');
% end
% disp('approx:');
% disp(max(n1));

r1 = zeros(16,16);
for i1 = 0:15
    r1(i1+1,:) = bitxor(0:15, i1*ones(1,16));
end
r1
    
    
    
    
    