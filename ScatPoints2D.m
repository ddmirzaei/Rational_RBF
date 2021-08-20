function [X,Xcov,PatchNearBound] = ScatPoints2D(a,b,h,hcov,PointType)
% X: interpolation points of size (N x d) 
% Xcov: covering centers of size (Nc x d)
% PatchNearBound: indices of patch centers which are close to the boundary 
% a,b: Omega = [a,b]^2
% h: spacing distance of points in X
% hcov: spacing distance of centers in Xcov
% PointType: the type of point distribution, here 'grid' or 'halton' 
switch (PointType)
  case ('grid')
    [xt,yt]=meshgrid(a:h:b,a:h:b); X=[xt(:) yt(:)];                      
  case('halton')
    N = round(((b-a)/h+1)^2); H = haltonset(2,'Skip',1e4,'Leap',1e2);
    X = net(H,N); X = (b-a)*X-a;
end
[xc,yc]=meshgrid(a:hcov:b,a:hcov:b);
xc1=xc(:); xc2=yc(:); Xcov = [xc1 xc2];
PatchNearBound = (xc1<a+hcov | xc1>b-hcov | xc2<a+hcov | xc2>b-hcov);
