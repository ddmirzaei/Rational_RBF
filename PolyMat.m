function P = PolyMat(X)
% X: points coordinates of size (N x d)
% P: Vandermonde matrix of size (N x Q) on monomials at X
% Note: This function works in 2D for polynomial orders <= 5, 
%   to generalize to 3D and higher orders just change the MultiIndex vector
global PolyOrder 
MultiIndex= [0 0;1 0;0 1;2 0;1 1;0 2;3 0;2 1;1 2;0 3;4 0;3 1;2 2;1 3;0 4]; 
d = size(X,2); 
Q = nchoosek(PolyOrder-1+d,d); % dimension of polynomial space in R^d
P = [];
for k = 1:Q
  P = [P prod(X.^MultiIndex(k,:),2)];
end