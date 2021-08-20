function K = KerMat(Xe,X)
% Xe: evaluation points of size (M x d)
% X:  trial points of size (N x d)
% K: kernel matrix of size (M x N) based on polyharmonic splines (PHS)
global RBFtype RBFpar
r = DistanceMat(Xe,X);
switch RBFtype
  case 'p'   % powers
    K = (-1).^ceil(RBFpar/2)*r.^RBFpar;                   
  case 'tps' % thin plate splines
    K = (-1)^(RBFpar/2+1)*(r+eps).^RBFpar.*log(r+eps); 
end