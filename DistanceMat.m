function DM = DistanceMat(Xe,X)
% Xe: evaluation points of size (M x d)
% X: trial points of size (N x d)
% DM: distance matrix of size (M x N)
[n,dim] = size(X); m = size(Xe,1); DM = 0;
for d = 1:dim
  DM = DM + (repmat(X(:,d)',m,1) - repmat(Xe(:,d),1,n)).^2;
end
DM = sqrt(DM);
