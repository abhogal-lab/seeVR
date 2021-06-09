function [p, err]=legendreP(l,m,x)

% Calculate the Associated Legendre polynomials. Speed is obtained by
% direct calculation of polynomial coefficients rather than recursion.
% Polynomial coefficients can increase in magnitude very quickly with
% polynomial degree, leading to decreased accuracy (estimated by err return
% parameter). If you need higher degrees for the polynomials, use
% recursion-based algorithms.

% If called with two return values, the second is a scalar, estimating the
% error, based on worst case scenario and the significant digits. Error is
% based on the largest polynomial coefficient and machine error, "eps"
% maxcf is used to store the largest polynomial coefficient If there are
% only two input arguments, the second argument is x, and m defaults to 0;
% then legendreP is the Legendre polynomial of degree m 
if(nargin==2), 
    x=m;
    clear m;
    m=0;
end;


%Some basic error checking of parameters

% Definition states that if |m|>l, the polynomial is 0; so set to 0 rather
% than return an error. Some algorithms depend on this behaviour. Note that
% the original definition requires 0<=m<=l
if(abs(m)>abs(l)), p=zeros(size(x)); err=0; return; end;

% Definition for l<0
if(l<0), l=-l-1; end;

% For m<0, polynomials are proportional to those with m>0, cfnm is the
% proportionality coefficient
cfnm=1;
if(m<0),
    m=-m;
    cfnm=(-1)^m*factorial(l-m)/factorial(l+m);
end;



% Calculate coef of maximum degree in x from the explicit analytical
% formula
cl=(-1)^m*cfnm*factorial(2*l)/((2^l)*factorial(l)*factorial(l-m));
maxcf=abs(cl);
%fprintf('Coef is %.16f\n',cl);
px=l-m;

% Power of x changes from one term to the next by 2. Also needed for
% sqrt(1-x^2).
x2=x.*x;


% Calculate efficiently P_l^m (x)/sqrt(1-x^2)^(m/2) - that is, only the
% polynomial part. At least one coefficient is guaranteed to exist - there
% is no null Legendre polynomial.
p=cl*ones(size(x));

for j=l-1:-1:0,
    % Check the exponent of x for current coefficient, px. If it is 0 or 1,
    % just exit the loop
    if(px<2), break; end;
    % If current exponent is >=2, there is a "next" coefficient; multiply p
    % by x2 and add it. Calculate the current coefficient
    cl=-(j+j+2-l-m)*(j+j+1-l-m)/(2*(j+j+1)*(l-j))*cl;
    
    if(maxcf<abs(cl)), maxcf=abs(cl); end;
    %fprintf('Coef is %.16f\n',cl);
    % ...and add to the polynomial
    p=p.*x2 + cl;
    % Decrease the exponent of x - this is the exponent of x corresponding
    % to the newly added coefficient
    px=px-2;
end;
% Estimate the error
err=maxcf*eps;
fprintf('Coef is %.16f, err %.16f\n',maxcf, err);

% Now we're done adding coefficients. However, if the exponent of x
% corresponding to the last added coefficient is 1 (polynomial is odd),
% multiply the polynomial by x 
if(px==1), p=p.*x; end;

% All that's left is to multiply the whole thing with sqrt(1-x^2)^(m/2). No
% further calculations are needed if m=0.
if(m==0), return; end;

x2=1-x2;
%First, multiply by the integer part of m/2
for j=1:floor(m/2), p=p.*x2; end;
%If m is odd, there is an additional factor sqrt(1-x^2)
if(m~=2*floor(m/2)), p=p.*sqrt(x2); end;

% Finally, the polynomials are not defined for |x|>1. If you do not need
% this behaviour, comment the following line
p(abs(x)>1)=NaN;

% Done.
