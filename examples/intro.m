%% Load the required modules
mrstModule add vortups coarsegrid mimetic incomp agglom upscaling
verbose = true;

%% Define fine-scale model and coarse grid

cellDims  = [50 50 1];
G         = cartGrid(cellDims, cellDims);
G         = computeGeometry(G);
clf, plotGrid(G);
load perm;
rock.perm = convertFrom(K(:), milli()*darcy());
rock.poro = repelem(0.2, numel(rock.perm))';
pv        = poreVolume(G, rock);
clf, 
plotGrid(G, 'EdgeColor', 'none');
plotCellData(G, log10(rock.perm(:,1)), 'EdgeColor', 'none'), axis equal tight off;

W     = addWell([], G, rock, 1,      ...
            'Type', 'rate' , 'Val', 100*meter^3/day(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1 0]);
W     = addWell(W, G, rock, G.cells.num,      ...
            'Type', 'bhp' , 'Val', 200*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0 1]);  
fluid = initSimpleFluid('mu' , [   1,   1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
iVor  = sphvormap(G, rock, W);
clf,
plotGrid(G, 'EdgeColor', 'none');
plotWell(G, W, 'radius', 2.5, 'color', 'w');
plotCellData(G, iVor, 'EdgeColor', 'none'), axis equal tight off;

% The coarse grid
coarsed = [10 10 1];
p0      = partitionUI(G, coarsed);
p1      = refineUniform(p0, G, iVor, 26, 'cartDims', [3 3 1]);
p1      = compressPartition(p1);
[bl,p1] = findConfinedBlocks(G, p1);
coarseG = generateCoarseGrid(G, p1);
coarseG = coarsenGeometry(coarseG);
pvC     = accumarray(coarseG.partition  ,pv);
 
outlineCoarseGrid(G, coarseG.partition, 'EdgeColor', 'k');
                         