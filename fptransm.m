function out = fptransm(f,param)
%fptransm(f,param) calculates the intensity transmittance of a Fabry-Perot cavity.
%Parameters: 
%   f: input frequency
%   param(1): peak transmision
%   param(2): free-spectral range
%   param(3): finesse

%here f is the frequency, not angular frequency
out = param(1)./(1+(2*param(3)/pi).^2*sin(pi*f/param(2)).^2);