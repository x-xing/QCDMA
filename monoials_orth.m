n=3;
phi = zeros(n);

phi(1,:)=(0:(n-1))*2*pi/n;
for ind = 2:n
    phi(ind,:)=circshift(phi(ind-1,:),[0 1]);
end
p1=exp(i*phi);
for ind1 = 1:n
    for ind2 = 1:n
        if abs(sum(p1(ind1,:).*p1(ind2,:)))>1e-8
            error('basis vector are not orthogonal!');
        end
    end
end
plot(p1(3,:))