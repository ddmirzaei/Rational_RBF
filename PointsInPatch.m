function Ind = PointsInPatch(X,Xcov)
% X: a set of points 
% Xcov: patch's centers
% Ind: indices of points in each patch
global PatchRadius PatchNearBound FacBound
T = KDTreeSearcher(X);   
Ind = cell(size(Xcov,1),1);
IndDom = rangesearch(T,Xcov(~PatchNearBound,:),PatchRadius);
IndBound = rangesearch(T,Xcov(PatchNearBound,:),FacBound*PatchRadius);
Ind(PatchNearBound) = IndBound;
Ind(~PatchNearBound) = IndDom;
