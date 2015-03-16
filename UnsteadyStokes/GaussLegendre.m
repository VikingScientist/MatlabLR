function [p,w] = GaussLegendre(np);

% This method takes in the number of points (np) associated with the
% GaussLegendre quadrature, and returns an np-dimensional array of
% the points (p) and weights (w).
%
% authors: Kay Hansen-Zahl
%          Tormod Bjøntegaard
% 
%          Fall 2002

A = zeros(np);
p = zeros(np,1);
w = zeros(np,1);

% This loop finds the A-matrix
A(1,2) = 1;
if np>2
  for i = 2:np-1
    A(i,i-1) = (i-1)/(2*i-1);
    A(i,i+1) = i/(2*i-1);
  end
end  
A(np,np-1) = (np-1)/(2*np-1);

% The array of the sorted eigenvalues/zeros
p=sort(eig(A));

% This loop finds the associated weights
for j=1:np
  w(j) = 2/((1-p((j))^2)*(LegendreDerivative(np,p(j)))^2);
end


