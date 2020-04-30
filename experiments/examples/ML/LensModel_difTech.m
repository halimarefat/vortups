clc;
clear;
close all;
%% Load the required modules
mrstModule add mrst-gui vortups coarsegrid mimetic incomp agglom upscaling 
verbose = true;
%% Cartesian fine Geometry is defined and initial solution is done
res_n = 1;
dims  = res_n * [80, 20, 1];
cells = dims(1,1)* dims(1,2);
G = cartGrid(dims, dims);
G = computeGeometry(G);
rock.perm = zeros(cells, 1);
rock.perm(:) = 1;
lens_dim1 = res_n * [48, 11, 1];
row = res_n * 6;
col = res_n * 10;
cell_str = row*dims(1,1)+ col;
cell_end = cell_str + lens_dim1(1,1);
lens_indx1 = zeros(lens_dim1(1,1)*lens_dim1(1,2), 1);
ii = 1;
for i = 1:lens_dim1(1,2)+1
    for j = 0:lens_dim1(1,1)
        lens_indx1(ii) = cell_str + j + (i-1)*dims(1,1);
        ii = ii+1;
    end 
end
rock.perm(lens_indx1(:)) = 0.001;

clf, 
plotGrid(G, 'EdgeColor', 'none');
plotCellData(G, log10(rock.perm(:,1)), 'EdgeColor', 'k'), axis equal tight off;
colorbar;
hold on

rock.poro = zeros(cells, 1);
rock.poro(:) = 0.2;

inj   = (1041:80:G.cells.num);
prod  = (80:80:600);
inj_s = [1,0;1,0;1,0;1,0;1,0;1,0;1,0];
prod_s= [0,1;0,1;0,1;0,1;0,1;0,1;0,1];
src   = addSource( [],  inj, 1000*meter^3/day(), 'sat', inj_s);
src   = addSource(src, prod,-1000*meter^3/day(), 'sat',prod_s);

MR    = 1; % mobility ratio     
fluid = initSimpleFluid('mu' , [   1,MR*1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);         
hT    = computeTrans(G, rock);
rSol  = initState(G, [], 0, 0);
S     = computeMimeticIP(G, rock);
rSol  = incompMimetic(rSol, G, S, fluid, 'src', src);
vor   = vorticitycalculator(G, rSol.flux); 
pv    = poreVolume(G, rock);

 vel = faceFlux2cellVelocity(G, rSol.flux);
 dim = [G.cartDims(1), G.cartDims(2)];
 u   = reshape(vel(:,1), dim);
 v   = reshape(vel(:,2), dim);
 X   = reshape(G.cells.centroids(:,1), dim);
 Y   = reshape(G.cells.centroids(:,2), dim);
quiver(X, Y, u, v, 'color', 'k');
text(G.cells.centroids(inj,1), G.cells.centroids(inj,2), int2str(G.cells.indexMap(inj,1)), 'color', 'r');
text(G.cells.centroids(prod,1), G.cells.centroids(prod,2), int2str(G.cells.indexMap(prod,1)), 'color', 'r');

plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
colorbar;
title('Initial pressure'), view(3)

cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
plotCellData(G, accumarray(cellNo, ...
    abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day))));
colorbar;
title('Initial flux intensity'), view(2)

plotCellData(G, abs(vor)); colorbar;
caxis([min(abs(vor)) max(abs(vor))])
title('Vorticity Map'), view(2)

iPer= log10(rock.perm(:,1));
iPer= iPer - min(iPer) + 1;

vel = sqrt(sum(vel .^ 2, 2));
iVel= log10(vel);
iVel= iVel - min(iVel) + 1;

iVor= abs(vor);  

p1  = partitionUI(G, [10, 4, 1]);

pPer= refineUniform(p1, G, iPer, 35, 'cartDims', [2 2 1]);
pPer= compressPartition(pPer);
[blksPer, pPer] = findConfinedBlocks(G, pPer);
CGPer = generateCoarseGrid(G, pPer);
CGPer = coarsenGeometry(CGPer);
pvCPer= accumarray(CGPer.partition  ,pv);
CGPer.cells.volumes = accumarray(CGPer.partition, G.cells.volumes);
CGPer.nodes.coords  = zeros(CGPer.cells.num, 3);
CGPer.faces.normals = zeros(CGPer.faces.num, 3);
[nsubCPer, subCPer] = subFaces(G, CGPer);
[sgnCPer, cfCPer] = signOfFineFacesOnCoarseFaces(G, CGPer, nsubCPer, subCPer);
rockCPer.perm = upscalePerm(G, CGPer, rock, 'Verbose', verbose);
rockCPer.poro = accumarray(CGPer.partition, pv)./CGPer.cells.volumes;
srcCPer       = convertSource2Coarse(CGPer, src);

pVel= refineUniform(p1, G, iVel, 35, 'cartDims', [2 2 1]);
pVel= compressPartition(pVel);
[blksVel, pVel] = findConfinedBlocks(G, pVel);
CGVel = generateCoarseGrid(G, pVel);
CGVel = coarsenGeometry(CGVel);
pvCVel= accumarray(CGVel.partition  ,pv);
CGVel.cells.volumes = accumarray(CGVel.partition, G.cells.volumes);
CGVel.nodes.coords  = zeros(CGVel.cells.num, 3);
CGVel.faces.normals = zeros(CGVel.faces.num, 3);
[nsubCVel, subCVel] = subFaces(G, CGVel);
[sgnCVel, cfCVel] = signOfFineFacesOnCoarseFaces(G, CGVel, nsubCVel, subCVel);
rockCVel.perm = upscalePerm(G, CGVel, rock, 'Verbose', verbose);
rockCVel.poro = accumarray(CGVel.partition, pv)./CGVel.cells.volumes;
srcCVel       = convertSource2Coarse(CGVel, src);

pVor= refineUniform(p1, G, iVor, 35, 'cartDims', [2 2 1]);
pVor= compressPartition(pVor);
[blksVor, pVor] = findConfinedBlocks(G, pVor);
CGVor = generateCoarseGrid(G, pVor);
CGVor = coarsenGeometry(CGVor);
pvCVor= accumarray(CGVor.partition  ,pv);
CGVor.cells.volumes = accumarray(CGVor.partition, G.cells.volumes);
CGVor.nodes.coords  = zeros(CGVor.cells.num, 3);
CGVor.faces.normals = zeros(CGVor.faces.num, 3);
[nsubCVor, subCVor] = subFaces(G, CGVor);
[sgnCVor, cfCVor] = signOfFineFacesOnCoarseFaces(G, CGVor, nsubCVor, subCVor);
rockCVor.perm = upscalePerm(G, CGVor, rock, 'Verbose', verbose);
rockCVor.poro = accumarray(CGVor.partition, pv)./CGVor.cells.volumes;
srcCVor       = convertSource2Coarse(CGVor, src);

close all;
figure;
plotCellData(G, rock.perm(:,1), 'EdgeColor', 'none'), axis equal tight off; hold on 
outlineCoarseGrid(G, CGVor.partition, 'Edgecolor', 'r', 'LineWidth', 0.8); 
figure;
plotCellData(G, log10(abs(vor)), 'EdgeColor', 'none'), axis equal tight off; hold on 
plotGrid(G, 'facecolor','none', 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5); 
outlineCoarseGrid(G, CGVor.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1); 
figure;
plotCellData(CGVor, rockCVor.perm(:,1), 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5), axis equal tight off; 
outlineCoarseGrid(G, CGVor.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1);


rcCPer.flux= accumarray(cfCPer, sgnCPer.*rSol.flux(subCPer), [CGPer.faces.num,1]);
rcCPer.s   = coarse_sat(rSol.s, CGPer.partition, pv, CGPer.cells.num);
rcPer      = deal(rSol);

rcCVel.flux= accumarray(cfCVel, sgnCVel.*rSol.flux(subCVel), [CGVel.faces.num,1]);
rcCVel.s   = coarse_sat(rSol.s, CGVel.partition, pv, CGVel.cells.num);
rcVel      = deal(rSol);

rcCVor.flux= accumarray(cfCVor, sgnCVor.*rSol.flux(subCVor), [CGVor.faces.num,1]);
rcCVor.s   = coarse_sat(rSol.s, CGVor.partition, pv, CGVor.cells.num);
rcVor      = deal(rSol);

close all;
step= 5;
T   = 3000;
dT  = T/step;
Tvec= 0:dT:T;
t   = 0; 
err = zeros(step+1, 3); err(2:end,:) = NaN;
p0  = get(gcf, 'OuterPosition');
clf, set(gcf, 'OuterPosition', [p0(1)-500, p0(2)-150, 1000, 600]);
mytitle =@(x) title(x, 'FontSize', 9, 'FontWeight', 'normal');
pos1 = 0.00;
pos2 = 0.15;
pos3 = 0.30;
pos4 = 0.75;
dis  = 0.1;
ind  = 0;
while t < T,
    ind  = ind + 1;
    
    rSol = implicitTransport(rSol, G, dT,  rock, fluid, 'src',  src);
    
    rcCPer  = implicitTransport(rcCPer, CGPer, dT, rockCPer, fluid, 'src', srcCPer);
	rcPer.s = rcCPer.s(CGPer.partition);
    
    rcCVel  = implicitTransport(rcCVel, CGVel, dT, rockCVel, fluid, 'src', srcCVel);
    rcVel.s = rcCVel.s(CGVel.partition);
    
    rcCVor  = implicitTransport(rcCVor, CGVor, dT, rockCVor, fluid, 'src', srcCVor);
    rcVor.s = rcCVor.s(CGVor.partition);
    
    clf
    colormap jet;      
      axes('position',[pos1 pos2 pos3 pos4])
      plotCellData(G, rSol.s, 'EdgeColor', 'none'), axis equal tight off
      view(90,90);
      mytitle(sprintf('Fine: %d cells', G.cells.num));
      
      axes('position',[pos1+dis pos2 pos3 pos4])
      plotCellData(G, rcPer.s, 'EdgeColor', 'none'), axis equal tight off
      outlineCoarseGrid(G, CGPer.partition, 'EdgeColor', 'w', 'EdgeAlpha', 0.3);
      view(90,90);
      mytitle(sprintf('Per: %d blocks', CGPer.cells.num));
      
      axes('position',[pos1+2*dis pos2 pos3 pos4])
      plotCellData(G, rcVel.s, 'EdgeColor', 'none'), axis equal tight off
      outlineCoarseGrid(G, CGVel.partition, 'EdgeColor', 'w', 'EdgeAlpha', 0.3);
      view(90,90);
      mytitle(sprintf('Vel: %d blocks', CGVel.cells.num));
      
      axes('position',[pos1+3*dis pos2 pos3 pos4])
      plotCellData(G, rcVor.s, 'EdgeColor', 'none'), axis equal tight off
      outlineCoarseGrid(G, CGVor.partition, 'EdgeColor', 'w', 'EdgeAlpha', 0.3);
      view(90,90);
      mytitle(sprintf('Vor: %d blocks', CGVor.cells.num));
      
      colorbar('southoutside', 'position',[pos1+0.08 pos2-0.05 4*dis-0.035 0.02])
      
      axes('position',[pos1+5.3*dis pos2+0.46 .4 .30])
      err(ind,:) = [sum(abs(rcPer.s - rSol.s).*pv), ...
                    sum(abs(rcVel.s - rSol.s).*pv), ...
                    sum(abs(rcVor.s - rSol.s).*pv)]./sum(rSol.s .* pv);
      plot(Tvec, err),
      legend('Per', 'Vel', 'Vor');
      set(gca,'XLim',[Tvec(1) Tvec(end)]), mytitle('Error')
    drawnow
    
    rSol = incompMimetic(rSol, G, S, fluid, 'src', src);
    
    rcPer   = incompMimetic(rcPer  , G, S, fluid, 'src', src);
    rcCPer.flux = accumarray(cfCPer, sgnCPer .* rcPer.flux(subCPer), [CGPer.faces.num, 1]);
    
    rcVel   = incompMimetic(rcVel  , G, S, fluid, 'src', src);
    rcCVel.flux = accumarray(cfCVel, sgnCVel .* rcVel.flux(subCVel), [CGVel.faces.num, 1]);
    
    rcVor   = incompMimetic(rcVor  , G, S, fluid, 'src', src);
    rcCVor.flux = accumarray(cfCVor, sgnCVor .* rcVor.flux(subCVor), [CGVor.faces.num, 1]);
    
    t = t + dT;
    if (t < T), continue, end
      
end
display('Done!')

