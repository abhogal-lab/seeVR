function [ L2 ] = LegendreN( n, x )
% This function returns Legendre polynomial of degree n.
%
% Input parameters:
% 	n - order of polynomial
%	x - vector of time axis
%
% Returned value:
%	L2 - Legendre polynomial of degree n

L1 = zeros( 1, length(x) ); 
L2 = ones( 1, length(x) );
for i = 1:n
    L0 = L1;
    L1 = L2;
    L2 = ( (2.*i-1) .* x .* L1 - (i-1) .* L0) ./ i;
end
end

