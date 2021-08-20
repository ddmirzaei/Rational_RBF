function sigma = RBF_PU(Xe, X, Xcov, fX, rho)
% This function implements the RBF-PU method with constant-generated weight
% Xe: evaluation points of size (Ne x d)
% X: interpolation (trial) points of size (N x d)
% Xcov: center of circular coverings of size (Nc x d)
% fX: f values at X of size (N x 1)
% rho: vector of patch radiuses
% sigma: output rational interpolant at Xe of size (Ne x 1)
[Nc,dim] = size(Xcov); Ne = size(Xe,1); sigma = zeros(Ne,1);
IndX = PointsInPatch(X,Xcov);   % indices of trial points in each patch
IndXe = PointsInPatch(Xe,Xcov); % indices of eval points in each patch
T = KDTreeSearcher(Xcov);       % indices of centers around each center  
IndXc = rangesearch(T,Xcov,(sqrt(dim)+0.1)*min(rho));   
for j = 1:Nc
  oj = Xcov(j,:);            % j-th patch center    
  Xj = X(IndX{j},:);         % points in the j-th patch
  inds0 = IndXe{j};          % indices of all eval points in the j-th patch 
  Xej0 = Xe(inds0,:);        % eval points in the j-th patch 
  D = DistanceMat(Xej0,Xcov(IndXc{j},:)); 
  [Dmin,Dind] = min(D,[],2); % find minimum distances
  % indices of eval points in which oj is their closest center
  inds = (Dind == 1);  
  indx = inds0(inds);        % going from local indices to global ones
  if ~isempty(inds)
    Xej = Xej0(inds,:);      % eval points those oj is their closest center
    fXj = fX(IndX{j});       % function values in the j-th patch
    sigmaj = Ratinal_RBF(Xej,Xj,oj,rho(j),fXj); % local interpolant
    sigma(indx) = sigmaj;    % putting local values in the global vector
  end
end
