function out = entropy(p_bias,T)
px = [p_bias;1-p_bias];
out = -sum(px.*log(px)/log(2).*([1;T]*ones(size(p_bias))),1);

