function Ixy = mutinfo_c(Px,Py,Pyx)
%mutinfo_c(Px,Py,Pyx) calculates the mutual infomation between random
%variable X and Y, given the conditional probability P(Y|X) and marginal
%P(X) and P(Y).

onesPx = ones(size(Px))';
%H(Y)
Hy = -sum(Py.*log(Py)/log(2));
%P(X,Y) = P(Y|X)*P(X)
Pxy = onesPx*Px.*Pyx;
%H(Y|X)
Hyx = sum(sum(Pxy.*log(onesPx*Py./Pxy)/log(2)));
Ixy = Hy - Hyx;
