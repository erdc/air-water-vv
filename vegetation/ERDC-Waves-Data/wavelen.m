function [L] = wavelen(d,T,n,g);
%
% function wavelength(d,T,n,g);
%
% n = # of iterations
%
% d = depth, this can be an array
% T = period
% n = number of iterations
% g = 32.2 or 9.81

Leck = (g*(T^2)*0.5/pi)*sqrt(tanh(4*pi*pi*d/(T*T*g)));  % 1984 SPM, p.2-7

L1 = Leck;

for k=1:n,
    L2 = (g*(T^2)*0.5/pi)*tanh(2*pi*d./L1);   % check if it's right
    L1 = L2;                                 % redo
end

L = L2;
k=2*pi/L;

ko = find(d<=0);
L(ko) = 0;