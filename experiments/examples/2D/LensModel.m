clc;
clear;
close all;
%% Load the required modules
mrstModule add mrst-gui vortups coarsegrid mimetic incomp agglom upscaling msrsb 
verbose = true;
%% Cartesian fine Geometry is defined and initial solution is done
res_n = 1;
dims  = res_n * [81, 30, 1];
cells = dims(1,1)* dims(1,2);
G = cartGrid(dims, dims);
G = computeGeometry(G);
rock.perm = zeros(cells, 1);
rock.perm(:) = 1;
lens_dim1 = res_n * [48,10, 1];
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
inj   =  1;
prod = G.cells.num;
W     = addWell(W, G, rock, inj,      ...
            'Type', 'rate' , 'Val', 500*meter^3/day(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W     = addWell(W, G, rock, prod,      ...
            'Type', 'bhp' , 'Val', 0*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]); 
plotWell(G, W);

MR    = 1; % mobility ratio      
fluid = initSimpleFluid('mu' , [   1,MR*1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);                     
hT    = computeTrans(G, rock);
Trans = hT2T(G, hT);
rSCart= initState(G, W, 100*barsa, [0, 1]);
SCart = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);
rSCart= incompMimetic(rSCart, G, SCart, fluid, 'Wells', W);
vor1   = vorticitycalculator(G, rSCart.flux);

pv = poreVolume(G, rock);

hT= computeTrans(G, rock);
Tr= hT2T(G, hT);
state0= initState(G,  W, 100*barsa, [0, 1]);
S = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);

                    
%% Pressure and Transport Equations are solved
totTime= 0.03*year; %1480657.84431263; %
N_step = 10;
dt     = totTime/N_step;

psolve = @(state) incompMimetic(state, G, S, fluid, 'Wells', W);
tsolve = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);

state  = psolve(state0); 
states = state;
mu_rf  = zeros(N_step, 2);
s_rf   = zeros(N_step, G.cells.num, 2);
kr_rf  = zeros(N_step, G.cells.num, 2);
lam_w_rf = zeros(N_step, G.cells.num);
lam_o_rf = zeros(N_step, G.cells.num);
F_w_rf = zeros(N_step, G.cells.num);
PVI_rf = zeros(N_step, 1);
for i = 1:N_step
    fprintf('Step %d of %d: ', i, N_step);
    
    fprintf('Solving pressure... ');
    state = psolve(states(end));
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states = [states; tsolve(state)];
    fprintf('Ok!\n');
    mu_rf(i,:)   = fluid.properties(states(end));    
    s_rf(i,:,:)  = fluid.saturation(states(end));
    for j = 1:G.cells.num
        kr_rf(i,j,:) = fluid.relperm(s_rf(i,j,:), states(end));
        lam_w_rf(i,j)= kr_rf(i,j,1) / mu_rf(i,1); % water mobility
        lam_o_rf(i,j)= kr_rf(i,j,2) / mu_rf(i,2); % oil mobility
    end
    F_w_rf(i,:) = lam_w_rf(i,:) ./ (lam_w_rf(i,:) + lam_o_rf(i,:));
    PVI_rf(i,1) = W(1).val(1) * (dt*i) / sum(pv.*rock.poro);
end
Fw_rf_prod = F_w_rf(:,prod);

p   = partitionUI(G, [8, 10]); %[6, 7, 1]

p   = compressPartition(p);
[blks, p] = findConfinedBlocks(G,p);
CG  = generateCoarseGrid(G, p);
CG  = coarsenGeometry(CG);
pvC = accumarray(CG.partition  ,pv);
close all;
figure;
plotCellData(G, rock.perm, 'EdgeColor', 'k'), axis equal tight off; hold on %rock.perm(:,1)
% colormap([1 1 1; parula]);
% plotGrid(G, 'facecolor','none', 'lineWidth', 0.25); view(-90,-90);
% outlineCoarseGrid(G, CGPEBI.partition, 'Color', 'k', 'lineWidth', 0.5); 
view(-90,-90);

%% Set up support regions required for MsRSB and move center points to wells
CG  = addCoarseCenterPoints(CG);
CG  = setCentersByWells(CG, W);
CG  = storeInteractionRegion(CG, 'ensureConnected', true);
% CGPEBI  = setupMexInteractionMapping(CGPEBI);

%% Set up basis functions
% (not now) By default, we use the C-accelerated version  
useCompiledBasis = true;
A = getIncomp1PhMatrix(G, hT);
% Update basis functions every now and then
updateBasis = true;

getBasis = @(A) getMultiscaleBasis(CG, A, 'type', 'MsRSB'); % ,'useMex', useCompiledBasis
basis0   = getBasis(A);

%% Solve multiscale 
% We solve the base case where only the multiscale solver is used.
basis = basis0;
W_ms  = W;
fn    = getSmootherFunction('type', 'ilu');

% without iterations
psolve = @(state, basis) incompMultiscale(state, CG, hT, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 0);
% with iterations
% psolve = @(state, basis) incompMultiscale(state, CGPEBI, hTPEBI, fluid, basis, 'wells', W_ms, ...
%     'getSmoother', fn, 'iterations', 5, 'useGMRES', true);

fprintf('Multiscale solution is started!\n');
states_ms = psolve(state0, basis);
mu_ms  = zeros(N_step, 2);
s_ms   = zeros(N_step, G.cells.num, 2);
kr_ms  = zeros(N_step, G.cells.num, 2);
lam_w_ms = zeros(N_step, G.cells.num);
lam_o_ms = zeros(N_step, G.cells.num);
F_w_ms = zeros(N_step, G.cells.num);
p0     = get(gcf, 'OuterPosition');
PVI_ms = zeros(N_step, 1);
err    = zeros(N_step, 1); err(2:end,:) = NaN;
clf, set(gcf, 'OuterPosition', [p0(1:2), 820, 640]);
mytitle =@(x) title(x, 'FontSize',10,'FontWeight','normal');
for i = 1:N_step
    state = states_ms(end);

    if updateBasis && mod(i, 10) == 0 && i > 1
        A = getIncomp1PhMatrix(G, hT, state, fluid);
        basis = getBasis(A);
    end
    
    fprintf('Step %d of %d: ', i, N_step);
    fprintf('Solving pressure... ');
    state = psolve(state, basis);
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_ms = [states_ms; tsolve(state)];
    fprintf('Ok!\n');
    
    mu_ms(i,:)   = fluid.properties(states_ms(end));    
    s_ms(i,:,:)  = fluid.saturation(states_ms(end));
    for j = 1:G.cells.num
        kr_ms(i,j,:) = fluid.relperm(s_ms(i,j,:), states_ms(end));
        lam_w_ms(i,j)= kr_ms(i,j,1) / mu_ms(i,1); % water mobility
        lam_o_ms(i,j)= kr_ms(i,j,2) / mu_ms(i,2); % oil mobility
    end
    F_w_ms(i,:) = lam_w_ms(i,:) ./ (lam_w_ms(i,:) + lam_o_ms(i,:));
    PVI_ms(i,1) = W_ms(1).val(1) * (dt*i) / sum(pv.*rock.poro);
    err(i,:)  = sum(abs(states_ms(i).s(:,1) - states(i).s(:,1)).*pv)./sum(states(i).s(:,1).*pv); 
%     % Plotting of saturations
%       clf
%       axes('position',[.04 .35 .2 .6])
%       plotCellData(G, states(i).s(:,1), 'EdgeColor', 'none'), axis equal tight off
%       mytitle(sprintf('Cart Fine: %d blocks', G.cells.num));           
% 
%       axes('position',[.28 .35 .2 .6])
%       plotCellData(G, states_ms(i).s(:,1), 'EdgeColor', 'none'), axis equal tight off
%       mytitle(sprintf('PEBI Fine: %d cells', G.cells.num));
%       
%       drawnow
end
Fw_ms_prod = F_w_ms(:,prod);

err      = zeros(N_step,1);
wtc      = zeros(N_step,2);
for i = 1:N_step
    err(i,1) = sum(abs(states_ms(i).s(:,1) - states(i).s(:,1)))./sum(states(i).s(:,1));
    wtc(i,:) = [states(i).s(prod,1), states_ms(i).s(prod,1)];   
end
names = {'Finescale', 'MsRSB'};

close all; plotToolbar(G, states);
axis equal tight off
daspect([1 1 0.2])
% view(85, 20);
% plotWell(G, W);
title(names{1});
colorbar('horiz')

figure; plotToolbar(G, states_ms);
axis equal tight off
% daspect([1 1 0.2])
% view(85, 20);
% plotWell(G, W);
title(names{2});
colorbar('horiz')

fprintf('Simulation is done!\n')
