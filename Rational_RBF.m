function sigma = Rational_RBF(Xe,X,xc,h,f)
% Xe: evaluation points in local domain
% X: interpolation points in local domain
% xc: center of the local domain
% h: scale number  
% f: vector of function values
% sigma: output rational interpolation at Xe
if norm(f,2) > 10*eps
  X = (X-xc)/h; Xe = (Xe-xc)/h;               % shifting and scaling                              
  P = PolyMat(X); [n,Q] = size(P);            % polynomial matrix
  Kmat = KerMat(X,X);                         % kernel matrix
  Z = null(P');                               % nullity of P^T
  ZKZinv = (Z'*Kmat*Z)\eye(size(Z,2));
  S = Z*ZKZinv*Z'; D = diag(f); 
  aa = 1/norm(f,2)^2; cc = aa; bb = 1; dd = 1;    
  A = aa*D*(S+eye(n))*D + bb*(S+eye(n)); 
  Binv = diag(1./sqrt(cc*f.^2+dd*ones(n,1)));
  A = Binv*A*Binv;
  [EigVec,EigVal] = eigs(A,1,'smallestabs');  % solving eigenvalue problem
  q = Binv*EigVec(:,end); p = D*q;            % q and p vectors
  KP = [Kmat P;P' zeros(Q)];                  % saddle point matrix
  Lmat = KP\[KerMat(Xe,X)';PolyMat(Xe)'];     % Lagrange functions at Xe
  p = Lmat(1:n,:)'*p; q = Lmat(1:n,:)'*q;     % numerator and denominator
  sigma = p./q;                               % rational interpolation
else % |f|_2 = 0
  sigma = zeros(size(Xe,1),1);     
end