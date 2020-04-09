%% Simulate Johansen 3D

clc
clear
close all 

if isnumeric(gcf)
    myZoom = @zoom;
else
    myZoom = @(varargin) [];
end

grdecl = fullfile('C:\Users\marefat\Desktop\Dena\DenaGRDplusPRPS.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G(1));
rock = grdecl2Rock(grdecl, G.cells.indexMap);
rock.poro = max(rock.poro, 1e-4);
rock.ntg  = max(rock.ntg, 1e-4); %this should not be ok!
rock.perm(:,2) = rock.perm(:,1);
rock.perm(:,3) = rock.perm(:,1);

figure;
plotCellData(G, log10(rock.perm(:, 1)));
axis equal tight off
daspect([1 1 0.2])
view(30, 10); myZoom(1.2);
colorbar
title('Dena Reservoir')
hold on

totTime = 50*year;
N_step = 100;
dt = totTime/N_step;

pv = poreVolume(G, rock);

wellIxI1 = [16, 16, 14, 14];
rate = 1.4e8*meter^3/day;
% W = verticalWell([], G, rock, wellIxI1(1), wellIxI1(2), wellIxI1(3):wellIxI1(4),...
%    'Type', 'rate', 'Val', rate, 'Radius', 0.1, 'comp_i', [1,0], ...
%    'name', 'I1', 'InnerProduct', 'ip_simple', 'sign', 1);
% wellIxI2 = [30, 30, 6, 6];
% W = verticalWell(W, G, rock, wellIxI2(1), wellIxI2(2), wellIxI2(3):wellIxI2(4),...
%    'Type', 'rate', 'Val', rate, 'Radius', 0.1, 'comp_i', [1,0], ...
%    'name', 'I2', 'InnerProduct', 'ip_simple', 'sign', 1);

wellIxP = [30, 30, 2, 2];
rateP = -1*1.4e8*meter^3/day;
W = verticalWell(W, G, rock, wellIxP(1), wellIxP(2), wellIxP(3):wellIxP(4),...
   'Type', 'rate', 'Val', rateP, 'Radius', 0.1, 'comp_i', [0,1], ...
   'name', 'P', 'InnerProduct', 'ip_simple', 'sign', -1);

plotWell(G, W);
hold off

T = getFaceTransmissibility(G, rock);
fluid = initSimpleFluid('mu', [1, 5]*centi*poise, 'n', [2, 2], 'rho', [0, 0]);

state0 = initResSol(G, 0, [0, 1]);
gravity reset off

psolve_TPFA  = @(state) incompTPFA(state, G, T, fluid, 'Wells', W, 'use_trans', true);
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

% cdims  = ceil(G.cartDims./[5, 5, 4]);
% padded = partitionUniformPadded(G, [cdims(1:2), 1]);
% uni = partitionUI(G, [1, 1, cdims(3)]);
% p1  = padded.*uni;
% p   = myrefineUniform(p1, G, iVor, 345, 'cartDims', [2 2 1]);
cdims  = ceil(G.cartDims./[6, 6, 2]);
padded = partitionUniformPadded(G, [cdims(1:2), 1]);
uni = partitionUI(G, [1, 1, cdims(3)]);
p   = padded.*uni;
p   = compressPartition(p);
fconn = ones(G.faces.num, 1);
fconn(G.faces.tag > 0) = 0;
p = mergeBlocksByConnections(G, p, fconn, 25);
p = processPartition(G, p);
p = compressPartition(p);

mp  = max(p);
CG  = generateCoarseGrid(G, p);%, cellPartitionToFacePartition(G,p));
CG  = coarsenGeometry(CG);
CG  = addCoarseCenterPoints(CG);
CG  = setCentersByWells(CG, W);
CG  = storeInteractionRegion(CG, 'ensureConnected', true);
CG  = setupMexInteractionMapping(CG);

useCompiledBasis = true;
A = getIncomp1PhMatrix(G, T);
updateBasis = true;
getBasis    = @(A) getMultiscaleBasis(CG, A, 'type', 'MsRSB');%, useCompiledBasis); % 'useMex',
basis0= getBasis(A);
basis = basis0;
W_ms  = W;
fn    = getSmootherFunction('type', 'ilu');

psolve_ms0it = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 0);
psolve_msItt = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 10, 'useGMRES', true);

state_ms0it  = psolve_ms0it(state0, basis);
states_ms0it = state_ms0it;

state_msItt  = psolve_msItt(state0, basis);
states_msItt = state_msItt;

err  = zeros(N_step, 2);
errp = zeros(N_step, 2);
% p0 = get(gcf, 'OuterPosition');
% clf, set(gcf, 'OuterPosition', [p0(1:2), 1500, 1040]);
mytitle =@(x) title(x, 'FontSize',10,'FontWeight','normal');
for i = 1:N_step
    fprintf('Step %d of %d: \n', i, N_step);
    
    fprintf('Solving FINEg pressure... ');
    state_TPFA = psolve_TPFA(states_TPFA(end));
    fprintf('Ok! ');
    fprintf('Solving FINEg transport... ');
    states_TPFA = [states_TPFA; solver(state_TPFA)];
    fprintf('Ok!\n');
    
    fprintf('Solving MS00t pressure... ');
    state_ms0it = psolve_ms0it(states_ms0it(end), basis);
    fprintf('Ok! ');
    fprintf('Solving MS00t transport... ');
    states_ms0it = [states_ms0it; solver(state_ms0it)];
    fprintf('Ok!\n');
    
    fprintf('Solving MS10t pressure... ');
    state_msItt = psolve_msItt(states_msItt(end), basis);
    fprintf('Ok! ');
    fprintf('Solving MS10t transport... ');
    states_msItt = [states_msItt; solver(state_msItt)];
    fprintf('Ok!\n');
    
    err(i,:) = [sum(abs(states_ms0it(i).s(:,1) - states_TPFA(i).s(:,1)).*pv), ...
                sum(abs(states_msItt(i).s(:,1) - states_TPFA(i).s(:,1)).*pv) ]./sum(states_TPFA(i).s(:,1) .* pv);
    
    errp(i,:)= [sum(abs(states_ms0it(i).pressure(:,1) - states_TPFA(i).pressure(:,1)).*pv), ...
                sum(abs(states_msItt(i).pressure(:,1) - states_TPFA(i).pressure(:,1)).*pv) ]./sum(states_TPFA(i).pressure(:,1) .* pv);
            
    clf
    plotCellData(G, states_msItt(i).s(:,1), 'EdgeColor', 'none'), axis equal tight off
    outlineCoarseGrid(G, CG.partition);
    mytitle(sprintf('Dena Res. Coarse: %d blocks', CG.cells.num));
    axis equal tight off
    daspect([1 1 0.2])
    view(30, 10); myZoom(1.2);
    colorbar
    plotWell(G, W);
    
    drawnow;
    
end



