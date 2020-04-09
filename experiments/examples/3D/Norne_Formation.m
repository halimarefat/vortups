%% Simulate Norne 3D

clc
clear
close all 

if ~(makeNorneSubsetAvailable() && makeNorneGRDECL()),
   error('Unable to obtain simulation model subset');
end

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G(1));
rock = grdecl2Rock(grdecl, G.cells.indexMap);

% if isnumeric(gcf)
%     myZoom = @zoom;
% else
%     myZoom = @(varargin) [];
% end
% 
% clf;
% plotGrid(G,'EdgeAlpha',.1,'FaceColor',[.9 .9 .7]); 
% view(70,70), zoom(2), axis tight off
% set(gca,'Zdir','normal'); camlight headlight; set(gca,'Zdir','reverse');
% 
% figure;
% plotCellData(G, log10(rock.perm(:, 1)));
% axis equal tight off
% daspect([1 1 0.2])
% view(85, 45); myZoom(1.2);
% colorbar
% title('Horizontal permeability (log10)')
% 
% figure;
% plotCellData(G, log10(rock.perm(:, 3)));
% axis equal tight off
% daspect([1 1 0.2])
% view(85, 45); myZoom(1.2);
% colorbar
% title('Vertical permeability (log10)')
% 
% 
% figure;
% plotCellData(G, rock.poro);
% axis equal tight off
% daspect([1 1 0.2])
% view(85, 45); myZoom(1.2);
% colorbar
% title('Porosity')

totTime = 100*year;
N_step = 100;
dt = totTime/N_step;

pv = poreVolume(G, rock);

wells = [13, 88,  -1; ...
         18, 87,  -1; ...
         36, 90,  -1; ...
         10, 15,  -1; ...
         24, 32,  1; ...
         8,  45,  1; ...
         16, 55,  1];

W = [];
[inum, pnum] = deal(1);
for i = 1:size(wells, 1);
    % Set well
    W = verticalWell(W, G, rock, wells(i, 1), wells(i, 2), [],...
                     'comp_i', [1, 0], 'type', 'bhp');
    if wells(i, 3) == 1
        % Producer
        W(i).val = -sum(pv)/(totTime*sum(wells(:, 3) == 1));
        W(i).type = 'rate';

        W(i).name = ['P', num2str(pnum)];
        W(i).sign = -1;
        pnum = pnum + 1;
    else
        % Injector
        W(i).val = 500*barsa;
        W(i).sign = 1;
        W(i).name = ['I', num2str(inum)];
        inum = inum + 1;
    end
end

% % Plot the grid, the wells and the perforated cells
% close all
% plotGrid(G, 'FaceColor', 'none', 'EdgeA', .2)
% plotWell(G, W)
% plotGrid(G, vertcat(W.cells), 'FaceColor', 'none', 'EdgeColor', 'b')
% axis equal tight off
% daspect([1 1 0.2])
% view(80, 65); myZoom(1.5);

T = getFaceTransmissibility(G, rock);
fluid = initSimpleFluid('mu', [1, 5]*centi*poise, 'n', [2, 2], 'rho', [0, 0]);

state0 = initResSol(G, 0, [0, 1]);
gravity reset off

psolve_TPFA  = @(state) incompTPFA(state, G, T, fluid, 'Wells', W, 'use_trans', true);
psolve_ms0it = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 0);
psolve_msItt = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 10, 'useGMRES', true);
solver       = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);

state_TPFA  = psolve_TPFA(state0);
vel    = faceFlux2cellVelocity(G, state_TPFA.flux); vel  = sqrt(sum(vel .^ 2, 2)); 
iVel   = log10(abs(vel));
iVel   = iVel - min(iVel) + 1;
states_TPFA = state_TPFA;

velUV  = faceFlux2cellVelocity(G, states_TPFA(1).flux);
vel    = sqrt(sum(velUV.^2,2));
visFl  = fluid.properties('mu');
gradmb = gradient(log(rock.perm(:,1)/visFl(1)));
c      = zeros(G.cells.num,3, 'gpuArray');
iVor   = zeros(G.cells.num,1, 'gpuArray');
for i = 1:G.cells.num
    a = [velUV(i,1), velUV(i,2), velUV(i,3)];
    b = [gradmb(i), gradmb(i), gradmb(i)];
    c(i,:)  = cross(a,b);
    iVor(i) = sqrt(c(i,1).^2 + c(i,2).^2 + c(i,3).^2);
end

% figure;
% plotCellData(G, log(iVor));
% axis equal tight off
% daspect([1 1 0.2])
% view(45, 45);
% colorbar
% title('Norne Res - Vorticity (ln)')
% hold on

cdims  = ceil(G.cartDims./[15, 10, 10]);
padded = partitionUniformPadded(G, [cdims(1:2), 1]);
uni = partitionUI(G, [1, 1, cdims(3)]);
p1  = padded.*uni;
% p   = refineUniformGPU(p1, G, iVor, 350, 'cartDims', [2 2 1]);
p   = myrefineUniform(p1, G, iVor, 350, 'cartDims', [2 2 1]);
    
% p1 = partitionUI(G, [11, 20, 11]);
% p  = refineUniform(p1, G, iVor, 300, 'cartDims', [2 2 1]);
% p  = processPartition(G, p);
% p  = compressPartition(p);
mp = max(p);
% 
CG = generateCoarseGrid(G, p);%, cellPartitionToFacePartition(G,p));
% plotCellData(G, log(iVor), 'EdgeColor', 'none')
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
CG = coarsenGeometry(CG);
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
% 
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
% 
CG = addCoarseCenterPoints(CG);
CG = setCentersByWells(CG, W);
CG = storeInteractionRegion(CG, 'ensureConnected', true);
CG = setupMexInteractionMapping(CG);

p0 = get(gcf, 'OuterPosition');
clf, set(gcf, 'OuterPosition', [p0(1:2), 1500, 1040]);
mytitle =@(x) title(x, 'FontSize',10,'FontWeight','normal');
for i = 1:N_step
    fprintf('Step %d of %d: ', i, N_step);
    
    fprintf('Solving pressure... ');
    state_TPFA = psolve_TPFA(states_TPFA(end));
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_TPFA = [states_TPFA; solver(state_TPFA)];
    fprintf('Ok!\n');
    
    clf
%     axes('position',[.01 .35 .2 .6])
    plotCellData(G, states_TPFA(i).s(:,1), 'EdgeColor', 'k'),     
    mytitle(sprintf('Fine: %d cells', G.cells.num));
    axis equal tight off    
    view(45, 45);
    
    
%     axes('position',[.16 .35 .2 .6])
%     plotCellData(G, rc.s, 'EdgeColor', 'none'), axis equal tight off
%     outlineCoarseGrid(G, CG.partition, 'Color', 'k', 'lineWidth', 1.5);
%     mytitle(sprintf('Coarse: %d blocks', CG.cells.num));

    drawnow;
    
end



