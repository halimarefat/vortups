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

iVor= abs(vor);  
p1  = partitionUI(G, [10, 4, 1]);
p   = refineUniform(p1, G, iVor, 35, 'cartDims', [2 2 1]);
p   = compressPartition(p);
[blks, p] = findConfinedBlocks(G, p);
CG  = generateCoarseGrid(G, p);
CG  = coarsenGeometry(CG);
pvC = accumarray(CG.partition  ,pv);

CG.cells.volumes = accumarray(CG.partition, G.cells.volumes);
CG.nodes.coords  = zeros(CG.cells.num, 3);
CG.faces.normals = zeros(CG.faces.num, 3);
[nsubC, subC] = subFaces(G, CG);
[sgnC, cfC] = signOfFineFacesOnCoarseFaces(G, CG, nsubC, subC);
rockC.perm = upscalePerm(G, CG, rock, 'Verbose', verbose);
rockC.poro = accumarray(CG.partition, pv)./CG.cells.volumes;
srcC  = convertSource2Coarse(CG, src);

close all;
figure;
plotCellData(G, rock.perm(:,1), 'EdgeColor', 'none'), axis equal tight off; hold on 
outlineCoarseGrid(G, CG.partition, 'Edgecolor', 'r', 'LineWidth', 0.8); 
figure;
plotCellData(G, log10(abs(vor)), 'EdgeColor', 'none'), axis equal tight off; hold on 
plotGrid(G, 'facecolor','none', 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5); 
outlineCoarseGrid(G, CG.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1); 
figure;
plotCellData(CG, rockC.perm(:,1), 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5), axis equal tight off; 
outlineCoarseGrid(G, CG.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1);
rcC.flux= accumarray(cfC, sgnC.*rSol.flux(subC), [CG.faces.num,1]);
rcC.s   = coarse_sat(rSol.s, CG.partition, pv, CG.cells.num);
rc      = deal(rSol);

T  = 5000;
dT = T/3;
t  = 0; 
p0 = get(gcf, 'OuterPosition');
clf, set(gcf, 'OuterPosition', [p0(1:2)-150, 820, 640]);
mytitle =@(x) title(x, 'FontSize',10,'FontWeight','normal');

while t < T,
    
    rSol = implicitTransport(rSol, G, dT,  rock, fluid, 'src',  src);
    rcC  = implicitTransport(rcC, CG, dT, rockC, fluid, 'src', srcC);
    rc.s = rcC.s(CG.partition);
    
    clf
      axes('position',[.04 .35 .2 .6])
      plotCellData(G, rSol.s, 'EdgeColor', 'none'), axis equal tight off
      view(90,90);
      mytitle(sprintf('Fine: %d cells', G.cells.num));
      axes('position',[.28 .35 .2 .6])
      plotCellData(G, rc.s, 'EdgeColor', 'none'), axis equal tight off
      outlineCoarseGrid(G, CG.partition, 'EdgeColor', 'w', 'EdgeAlpha', 0.3);
      view(90,90);
      mytitle(sprintf('Coarse: %d blocks', CG.cells.num));
     
    drawnow
    
    rSol = incompMimetic(rSol, G, S, fluid, 'src', src);
    rc   = incompMimetic(rc  , G, S, fluid, 'src', src);
    rcC.flux = accumarray(cfC, sgnC .* rc.flux(subC), [CG.faces.num, 1]);
    
    t = t + dT;
    if (t < T), continue, end
      
end
display('Done!')

