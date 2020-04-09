%% Simulate the SPE10 Model2 3D

clc
clear, close all hidden

spe10_data  = fullfile(fileparts(mfilename('fullpath')), ...
                       '..', 'spe10_rock.mat');
if ~exist(spe10_data, 'file'),
   if ~make_spe10_data,
      error(['Failed to establish on-disk representation of ', ...
             'SPE10 rock data']);
   end
end

%%
layers = 1:85;

cartDims = [60, 220, numel(layers)];
rock     = SPE10_rock(layers);

rock.perm = convertFrom(rock.perm, milli*darcy);
% An isotropic media
rock.perm(:,2) = rock.perm(:,1);
rock.perm(:,3) = rock.perm(:,1);

is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));

physDims = cartDims .* [20, 10, 2]*ft;

G = computeGeometry(cartGrid(cartDims, physDims));
S = computeMimeticIP(G, rock);

figure;
plotCellData(G, log10(rock.perm(:, 1)));
axis equal tight off
daspect([1 1 0.2])
view(45, 45);
colorbar
title('SPE10 Model 2 - permeability (log10)')
hold on

fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

well_ip = 'ip_simple'; %'ip_tpf';
% W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
%                  'InnerProduct', well_ip, ...
%                  'Val', 4000*psia, 'Radius', 0.125*meter, ...
%                  'Name', 'P1', 'Comp_i', [0, 0]);
% 
% W = verticalWell(W , G, rock, 60,   1, [], 'Type', 'bhp', ...
%                  'InnerProduct', well_ip, ...
%                  'Val', 4000*psia, 'Radius', 0.125*meter, ...
%                  'Name', 'P2', 'Comp_i', [0, 0]);

% W = verticalWell(W , G, rock, 60, 220, [], 'Type', 'bhp', ...
%                  'InnerProduct', well_ip, ...
%                  'Val', 4000*psia, 'Radius', 0.125*meter, ...
%                  'Name', 'P3', 'Comp_i', [0, 0]);

% W = verticalWell(W , G, rock,  1, 220, [], 'Type', 'bhp', ...
%                  'InnerProduct', well_ip, ...
%                  'Val', 4000*psia, 'Radius', 0.125*meter, ...
%                  'Name', 'P4', 'Comp_i', [0, 0]);
% 
% W = verticalWell(W , G, rock, 30, 110, [], 'Type', 'rate',   ...
%                  'InnerProduct', well_ip, ...
%                  'Val', 5000*stb/day, 'Radius', 0.125*meter, ...
%                  'Name', 'I1', 'Comp_i', [1, 0]);

W = verticalWell([] , G, rock, 1, 1, [], 'Type', 'rate',   ...
                 'InnerProduct', well_ip, ...
                 'Val', 5000*stb/day, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);
W = verticalWell(W , G, rock, 60, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 4000*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P2', 'Comp_i', [0, 0]);
             
plotWell(G, W);
hold off

totTime = 100*year;
N_step  = 100;
dt      = totTime/N_step;

pv = poreVolume(G, rock);
T  = getFaceTransmissibility(G, rock);

state0 = initResSol(G, 0, [0, 1]);
psolve = @(state) incompTPFA(state, G, T, fluid, 'Wells', W, 'use_trans', true);
solver = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);
state  = psolve(state0);
states = state;

% velUV = faceFlux2cellVelocity(G, states(1).flux);
% vel   = sqrt(sum(velUV.^2,2));
% visFl = fluid.properties('mu');
% gradmb= gradient(log(rock.perm(:,1)/visFl(1)));
% c     = zeros(G.cells.num,3);
% iVor  = zeros(G.cells.num,1);
% for i = 1:G.cells.num
%     a = [velUV(i,1), velUV(i,2), velUV(i,3)];
%     b = [gradmb(i), gradmb(i), gradmb(i)];
%     c(i,:)  = cross(a,b);
%     iVor(i) = sqrt(c(i,1).^2 + c(i,2).^2 + c(i,3).^2);
% end 
% iVel= log10(abs(vel));
% iVel= iVel - min(iVel) + 1;

% figure;
% plotCellData(G, log(iVor));
% axis equal tight off
% daspect([1 1 0.2])
% view(45, 45);
% colorbar
% title('SPE10 Model 2 - Vorticity (ln)')
% hold on
% 
% p1 = partitionUI(G, [6, 11, 5]);
% p  = refineUniform(p1, G, iVor, 5000, 'cartDims', [2 2 1]);
% p  = processPartition(G, p);
% p  = compressPartition(p);
% mp = max(p);
% 
% CG = generateCoarseGrid(G, p, cellPartitionToFacePartition(G,p));
% plotCellData(G, mod(p, 13), 'EdgeColor', 'none')
% plotGrid(CG, 'facec', 'none', 'edgec', 'k', 'linewidth', 2)
% hold off
% 
% figure;
% plotCellData(CG,(1:mp)');
% plotFaces(CG,1:CG.faces.num,'FaceColor' , 'none' , 'LineWidth' ,2);
% view(3); axis off
% colormap(.5*(jet(128)+ones(128,3)));
% hold on
% 
% CG = coarsenGeometry(CG);
% plotPts = @(pts, varargin) plot3(pts(:,1), pts(:,2), pts(:,3), varargin{:});
% hold on
% h=plotPts(CG.faces.centroids, 'k*');
% hold off;
% 
% i = false(mp,1); i(4) = true;
% cla, hold
% cg_cent = CG.cells.centroids(i,:);   plotPts(cg_cent, 'ok','MarkerFaceColor',[.5 .5 .5]);
% g_cent = G.cells.centroids(p==4,:);  plotPts(g_cent, '.');
% plotCellData(CG,(1:mp)',~i);
% plotFaces(CG,gridCellFaces(CG,find(~i)),'FaceColor','none','LineWidth',2);
% view(-35,15);
% legend({'Coarse centroids', 'Fine Centroids'}, 'Location', 'SouthOutside')

% face  = 66;
% sub   = CG.faces.connPos(face):CG.faces.connPos(face+1)-1;
% ff    = CG.faces.fconn(sub);
% neigh = CG.faces.neighbors(face,:);
% 
% figure;
% show = false(1,CG.faces.num);
% show(boundaryFaces(CG)) = true;
% show(boundaryFaces(CG,neigh)) = false;
% plotFaces(CG, show,'FaceColor',[1 1 .7]);
% plotFaces(CG,boundaryFaces(CG,neigh),'FaceColor','none','LineWidth', 2);
% plotFaces(G, ff, 'FaceColor', 'g')
% plotGrid(G, p == neigh(1), 'FaceColor', 'none', 'EdgeColor', 'r')
% plotGrid(G, p == neigh(2), 'FaceColor', 'none', 'EdgeColor', 'b')
% view(3); axis off
% 
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% 
% figure; 
% plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
% plotGrid(CG, 'facec', 'none', 'edgec', 'k', 'linewidth', 2);
% % plotCellData(G, log(iVor), j == round(G.cartDims(2)/2))
% plotCellData(G, log(iVor), abs(log(iVor))>21)
% axis equal tight off
% plotWell(G, W);
% view(45,50)

% CG = addCoarseCenterPoints(CG);
% CG = setCentersByWells(CG, W);
% CG = storeInteractionRegion(CG, 'ensureConnected', true);
% CG = setupMexInteractionMapping(CG);
% 
figure;
for i = 1:N_step
    fprintf('Step %d of %d: ', i, N_step);
    
    fprintf('Solving pressure... ');
    state = psolve(states(end));
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states = [states; solver(state)];
    fprintf('Ok!\n');
    
    plotCellData(G, states(i).s(:,1), 'EdgeColor', 'k'),     
    axis equal tight off
    daspect([1 1 0.2])
    view(45, 45);
    drawnow;
    
end



