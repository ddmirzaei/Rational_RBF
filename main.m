global RBFpar RBFtype PolyOrder 
global PatchNearBound FacBound PatchRadius 
RBFtype = 'tps'; RBFpar = 4;                  % polyharmonic RBF
PolyOrder = floor(RBFpar/2)+1;                % polynomial order 
Fexact = @(x,y) tan(9*(y-x)+1)/(tan(9)+1);    % function to be interpolated
he = 1/190; [xe,ye] = meshgrid(0:he:1,0:he:1);% mesh norm of eval points
Xe = [xe(:) ye(:)];                           % evaluation points
f_exact = Fexact(Xe(:,1),Xe(:,2));            % exact solution
FacBound = 1.5;         % incresing factor for patches near the boundary
h=0.05; niter = 8;
for k=1:niter
  hcov = 4*h;                            % spacing between patch centers
  [X,Xcov,PatchNearBound] = ScatPoints2D(0,1,h,hcov,'halton'); 
  N(k) = size(X,1);                      % number of interpolation points 
  C_ovlp = 1;                            % overlapping constant 
  PatchRadius = C_ovlp*hcov;             % basic patch radius
  rho = PatchRadius*ones(size(Xcov,1),1);% radius of internal patches
  rho(PatchNearBound)=FacBound*PatchRadius;% radius of near bondary patches
  fX = Fexact(X(:,1),X(:,2));            % f values at interpolation points
  sigma = RBF_PU(Xe,X,Xcov,fX,rho);      % rational RBF-PU interpolant    
  err(k,1) = norm(sigma-f_exact,2)/norm(f_exact,2) % interpoltion errors
  h = h/sqrt(2);                                   % point refinement
end
figure; 
loglog(N.^0.5,err,'-<b','MarkerSize',6,'MarkerFaceColor','b','LineWidth',1)
xlabel('$\sqrt{N}$','interpreter','latex'); xlim([15 500])
ylabel('$\|f-\sigma\|_{2,X_e}/\|f\|_{2,X_e}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[300 300 350 400]); set(gca,'XTick',[20 50 100 200 400])
