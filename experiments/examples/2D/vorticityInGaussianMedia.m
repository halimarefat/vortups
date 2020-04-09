clc;
clear;
close all;

res_n = 1;
nn= 80;
dims  = res_n * [nn, nn];
cells = dims(1,1)* dims(1,2);
GC = cartGrid(dims, dims);
GC = computeGeometry(GC);
permi  = abs(gaussianField(dims));
rockCa.perm = zeros(cells, 1);
rockCa.perm(:) = 1; %reshape(permi, [cells, 1]);
% G1C = cartGrid([20, 20]);
% G1C = computeGeometry(G1C);
% figure; 
% plotGrid(G1C);
% plotCellData(G1C,per,'EdgeColor','none'); axis tight off; axis equal tight
% 
% for i = 1:60
%     wIndx(i) = 40+(i-1)*80;
%     rockCa.perm(40+(i-1)*80) = 0.000000000000001; 
%     rockCa.perm(41+(i-1)*80) = 0.000000000000001; 
%     rockCa.perm(42+(i-1)*80) = 0.000000000000001; 
% end

% for i = 1:nn*nn/2
%      rockCa.perm(i) = 0.1; 
% end

% for i = 1:nn/2
%      rockCa.perm(40+(i-1)*nn) = 0.1; 
% end

for i = 1:nn*nn/2
     rockCa.perm(i) = 0.1; 
end

% for j = 0:nn-1
%     for i = (j*nn+1):(j*nn+41)
%          rockCa.perm(i) = 0.1; 
%     end
% end

rockCa.poro = zeros(cells, 1);
rockCa.poro(:) = 0.2;
%% Initialize model problem and compute initial solution
gravity off
fluid = initSimpleFluid('mu' , [   1,  .2]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
srcCa = addSource([], [1 GC.cells.num], [10000 -10000]./day(), 'sat', [1, 0]);

%% Compute initial fine-grid solution
hTC   = computeTrans(GC, rockCa);
TransC= hT2T(GC, hTC);
rSCart= initState(GC, [], 0, [0, 1]);
SCart = computeMimeticIP(GC, rockCa);
rSCart= incompMimetic(rSCart, GC, SCart, fluid, 'src', srcCa);
vor1  = vorticitycalculator(GC, rSCart.flux);
velxy = faceFlux2cellVelocity(GC, rSCart.flux);
vel1  = faceFlux2cellVelocity(GC, rSCart.flux);
vel1  = sqrt(sum(vel1.^2,2));
grdk  = gradient(log(rockCa.perm(:,1)));
lnK   = reshape(log(rockCa.perm(:,1)), [nn,nn]);
[gkx,gky]  = gradient(lnK, 1, 1);
xX = 1:1:nn;
yY = 1:1:nn;
contour(xX,yY,lnK)
hold on
quiver(xX,yY,gkx,gky)
hold off
view(-90,90);
% tempp = reshape(grdk, [nn,nn]);
vor2 = zeros(cells,1);
for i = 1:cells
    a = [velxy(i,1), velxy(i,2), 0];
    b = [gkx(i), gky(i), 0];
    c = cross(a,b);
    vor2(i) = c(3);
end 

figure; 
plotCellData(GC,log(rockCa.perm(:,1)/(milli*darcy)),'EdgeColor','none'); axis tight off;
plotGrid(GC, 'facecolor','none'); axis equal tight
title('CartGrid - Perm');

figure; 
plotCellData(GC,(abs(vor1)),'EdgeColor','none');  axis tight off; axis equal tight;
plotGrid(GC, 'facecolor','none'); axis equal tight
title('CartGrid - vor1');

figure; 
plotCellData(GC,(abs(vor2)),'EdgeColor','none');  axis tight off; axis equal tight;
plotGrid(GC, 'facecolor','none'); axis equal tight
title('CartGrid - vor2');
