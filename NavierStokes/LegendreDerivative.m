function Ld = LegendreDerivative(n,x);

% This function returns the value of the n'th derivative of the 
% Legendre polynomial evaluated in the point, x.
%
% authors: Kay Hansen-Zahl
%          Tormod Bjøntegaard
% 
%          Fall 2002

Ln = zeros(n+1,1);

Ln(1) = 1;
Ln(2) = x;

Ld(1) = 0;
Ld(2) = 1;

% Have to treat the endpoints separatly
if abs(x)==1
  Ld = x^(n-1)*(1/2)*n*(n+1);
else
  for i = 1:n-1
    Ln(i+2) = (2*i+1)/(i+1)*x*Ln(i+1) - i/(i+1)*Ln(i);
  end
  Ld = n/(1-x^2)*Ln(n) - n*x/(1-x^2)*Ln(n+1);
end
  


