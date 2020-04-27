clc;
clear;
close all;
%% Load the required modules
mrstModule add mrst-gui vortups coarsegrid mimetic incomp agglom upscaling msrsb steady-state
verbose = true;
%% Cartesian fine Geometry is defined and initial solution is done
res_n = 1;
dims  = res_n * [80, 20, 1];
cells = dims(1,1)* dims(1,2);
G = cartGrid(dims, dims);
G = computeGeometry(G);
rock.perm = zeros(cells, 1);
rock.perm(:) = 1;
lens_dim1 = res_n * [40, 10, 1];
row = res_n * 7;
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
hold on

rock.poro = zeros(cells, 1);
rock.poro(:) = 0.2;

W     = [];
% inj   = (1:80:1000);
% prod  = (1440:80:G.cells.num);
inj   = (561:80:G.cells.num);
prod  = (80:80:800);
for i = 1:numel(inj)
W     = addWell(W, G, rock, inj(i),      ...
            'Type', 'rate' , 'Val', 7.7*meter^3/day(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
end
for i = 1:numel(prod)
W     = addWell(W, G, rock, prod(i),      ...
            'Type', 'bhp' , 'Val', 0*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]); 
end
plotWell(G, W);
title('Permeability Map - Log'), view(3)

MR    = 1; % mobility ratio      
fluid = initSimpleFluid('mu' , [   1,MR*1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);                     
hT    = computeTrans(G, rock);
Trans = hT2T(G, hT);
rSol  = initState(G, W, 100*barsa, [0, 1]);
S     = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);
gravity off
rSol  = incompMimetic(rSol, G, S, fluid, 'Wells', W);
vor   = sphvormap(G, rock, W);
pv    = poreVolume(G, rock);

plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
title('Initial pressure'), view(3)

cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
plotCellData(G, accumarray(cellNo, ...
    abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day))));
title('Initial flux intensity'), view(3)

plotCellData(G, vor);
title('Vorticity Map'), view(2)

iVor= (vor); 
% iVor = iVor - min(iVor) + 1;  
p1  = partitionUI(G, [10, 4, 1]);
p   = refineUniform(p1, G, iVor, 41, 'cartDims', [2 2 1]);
p   = compressPartition(p);
[blks, p] = findConfinedBlocks(G, p);
CG  = generateCoarseGrid(G, p);
CG  = coarsenGeometry(CG);
pvC = accumarray(CG.partition  ,pv);
rockC.perm = upscalePerm(G, CG, rock, 'Verbose', verbose);

close all;
figure;
plotCellData(G, rock.perm(:,1), 'EdgeColor', 'none'), axis equal tight off; hold on 
plotGrid(G, 'facecolor','none', 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5); 
outlineCoarseGrid(G, CG.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1); 
figure;
plotCellData(G, log10(abs(vor)), 'EdgeColor', 'none'), axis equal tight off; hold on 
plotGrid(G, 'facecolor','none', 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5); 
outlineCoarseGrid(G, CG.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1); 
figure;
plotCellData(CG, rockC.perm(:,1), 'EdgeColor', [0.4,0.4,0.4], 'LineWidth', 0.5), axis equal tight off; 
outlineCoarseGrid(G, CG.partition, 'Edgecolor', [0.2,0.2,0.2], 'LineWidth', 1.1);

T      = 30;
dT     = T/3;
dTplot = 3;  % plot only every 100th day
N      = fix(T/dTplot);

rISol = rSol;
t  = 0; plotNo = 1; hi = 'Implicit: '; he = 'Explicit: ';
e = []; pi = []; pe = [];
while t < T,
   rSol  = explicitTransport(rSol , G, dT, rock, fluid, 'wells', W);
   rISol = implicitTransport(rISol, G, dT, rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   s = [rSol.s(:,1); rISol.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   rSol  = incompTPFA(rSol , G, hT, fluid, 'wells', W);
   rISol = incompTPFA(rISol, G, hT, fluid, 'wells', W);

   % Measure water saturation in production cells in saturation
   e = [e; sum(abs(rSol.s(:,1) - rISol.s(:,1)).*pv)/sum(pv)]; %#ok
   pe = [pe; rSol.s(W(2).cells,1)' ];                 %#ok
   pi = [pi; rISol.s(W(2).cells,1)'];                 %#ok

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t <T), continue, end

   % Plot saturation
   heading = [num2str(convertTo(t,day)),  ' days'];
   r = 0.01;
   subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
   plotCellData(G, rSol.s(:,1));
   view(60,50), axis equal off, title([he heading])

   subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
   plotCellData(G, rISol.s(:,1));
   view(60,50), axis equal off, title([hi heading])
   drawnow

   plotNo = plotNo+1;
end

%%
% As we clearly can see from the plots in the figure, the implicit scheme
% has much more numerical diffusion than the explicit scheme early in the
% simulation, but as the time increase, the difference is smaller. To
% verify this, we can plot the error or the breakthrough curves
%
n = size(pe,1);
pargs = {'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]};
subplot(2,1,1),
   plot(1:n,e*100,'-o', pargs{:}),
   title('Percentage saturation discrepancy')
subplot(2,1,2),
   plot(1:n,pe(:,1),'-o',1:n,pi(:,1),'-s',pargs{:})
   legend('Explicit','Implicit','Location','NorthWest');
   title('Water breakthrough at heel'); axis tight