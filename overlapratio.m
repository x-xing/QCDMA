function out = overlapratio(param, d1,d2,f)
% Overlap Ratio function overlapratio(d1,d2,filter) calculates the ratio of
% overlaps d1*filter(param) and d2*filter(param).
param(1)=0;
out = sum(d1.*Ldistr(f, param))/sum(d2.*Ldistr(f, param));