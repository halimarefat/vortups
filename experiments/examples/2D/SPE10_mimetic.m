clc;
clear;
close all;

%% Load the required modules
mrstModule add mrst-gui vortups coarsegrid mimetic incomp spe10 steady-state agglom upscaling
%% Cartesian fine Geometry is defined and initial solution is done
refined      = true;
[G, ~, rock] = getSPE10setup(56);
rock.poro    = max(rock.poro, 1e-4);

figure;
plotCellData(G, log10(rock.perm(:,1)), 'EdgeColor', 'none'), axis equal tight off;
view(90,90);

% adding well information and fluid properties
W     = [];
inj   = 1;
prodC = G.cells.num;
W     = addWell(W, G, rock, inj,      ...
            'Type', 'rate' , 'Val', 100*meter^3/day(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W     = addWell(W, G, rock, prodC,      ...
            'Type', 'bhp' , 'Val', 200*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]);  
fluid = initSimpleFluid('mu' , [   1,   1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

% Setting initial state of the reservoir                   
hT    = computeTrans(G, rock);
Trans = hT2T(G, hT);
S     = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);
pv    = poreVolume(G, rock);
state0= initState(G,  W, 100*barsa, [0, 1]);
%% Pressure and Transport Equations are solved for fine grid
totTime= 1.5*year;
N_step = 100;
dt     = totTime/N_step;

psolve = @(state) incompMimetic(state, G, S, fluid, 'Wells', W);
tsolve = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);

state    = psolve(state0); 
states   = state;
mu_rf    = zeros(N_step, 2);
s_rf     = zeros(N_step, G.cells.num, 2);
kr_rf    = zeros(N_step, G.cells.num, 2);
lam_w_rf = zeros(N_step, G.cells.num);
lam_o_rf = zeros(N_step, G.cells.num);
F_w_rf   = zeros(N_step, G.cells.num);
PVI_rf   = zeros(N_step, 1);
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
Fw_rf_prod = F_w_rf(:,prodC);

%% Cart Geometry Coarsening and Dual Mesh 
iVor = sphvormap(G, rock, W);

coa_dim = [3, 10, 1];
if refined
    p1 = partitionUI(G, coa_dim);
    p  = refineUniform(p1, G, iVor, 448, 'cartDims', [3 3 1]);
    p  = compressPartition(p);
else
    p  = partitionUI(G, coa_dim);
end
[blksCart, p] = findConfinedBlocks(G,p);
CG  = generateCoarseGrid(G, p);
CG  = coarsenGeometry(CG);
pvC = accumarray(CG.partition  ,pv);

CG.cells.volumes = accumarray(CG.partition, G.cells.volumes);
rockC      = convertRock2Coarse(G, CG, rock);
rockC.perm = upscalePerm(G, CG, rock, 'T', hT, 'S', S);
ChT        = computeTrans(CG, rockC);
CTrans1    = hT2T(CG, ChT);
[nsubC, subC] = subFaces(G, CG);
[sgnC, cfC  ] = signOfFineFacesOnCoarseFaces(G, CG, nsubC, subC);
[~,~,WC     ] = upscaleTrans(CG, hT, 'match_method', 'lsq_flux', ...
                             'bc_method', 'wells', 'wells', W); 
CS            = computeMimeticIP(CG, rockC, 'InnerProduct', 'ip_quasirt', ...
                                'FaceTrans', ...
                                [CG.cells.faces(:,1), ChT]);
figure;
plotCellData(G, log10(rock.perm(:,1)), 'Edgecolor', 'none'), axis equal tight off;
outlineCoarseGrid(G, CG.partition, 'Edgecolor', 'k', 'LineWidth', 0.5);
view(90,90);

%% Solve incompressible mimetic 
psolve = @(state) incompMimetic(state, CG, CS, fluid, 'Wells', WC);
tsolve = @(state) implicitTransport(state, CG, dt, rockC, fluid, 'wells', WC);

fprintf('Coarse solution is started!\n');
state0_cm= initState(CG,  WC, 100*barsa, [0, 1]);
states_cm= psolve(state0_cm);
mu_ms    = zeros(N_step, 2);
s_ms     = zeros(N_step, CG.cells.num, 2);
kr_ms    = zeros(N_step, CG.cells.num, 2);
lam_w_ms = zeros(N_step, CG.cells.num);
lam_o_ms = zeros(N_step, CG.cells.num);
F_w_ms   = zeros(N_step, CG.cells.num);
p_L2     = zeros(N_step, 1);
p_Linf   = zeros(N_step, 1);
f_L2     = zeros(N_step, 1);
f_Linf   = zeros(N_step, 1);
PVI_ms   = zeros(N_step, 1);
pres_tmp = zeros(G.cells.num, 1);
s_tmp    = zeros(G.cells.num, 1);
s_L1     = zeros(N_step, 1);
s_L2     = zeros(N_step, 1);
s_Linf   = zeros(N_step, 1);
for i = 1:N_step
    state = states_cm(end);
    fprintf('Step %d of %d: ', i, N_step);
    fprintf('Solving pressure... ');
    state = psolve(state);
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_cm = [states_cm; tsolve(state)];
    pres_tmp(:) = states_cm(i).pressure(p(:));
    s_tmp(:)  = states_cm(i).s(p(:), 1);
    p_L2(i)   = sqrt(sum(((pres_tmp - states(i).pressure).^2)) ./ sum((states(i).pressure.^2)));
    p_Linf(i) = max(abs(states(i).pressure - pres_tmp)) ./ max(abs(states(i).pressure));
    s_L1(i)   = sum((abs(s_tmp - states(i).s(:,1)).*pv)) ./ sum((states(i).s(:,1).*pv));
    s_L2(i)   = sqrt(sum(((s_tmp - states(i).s(:,1)).^2)) ./ sum((states(i).s(:,1).^2)));
    s_Linf(i) = max(abs(states(i).s(:,1) - s_tmp)) ./ max(abs(states(i).s(:,1)));
%     f_L2(i)   = sqrt(sum(((states_cm(i).flux - states(i).flux).^2)) ./ sum((states(i).flux.^2)));
%     f_Linf(i) = max(abs(states(i).flux - states_cm(i).flux)) ./ max(abs(states(i).flux));
    fprintf('Ok!\n');
%     mu_ms(i,:)   = fluid.properties(states_cm(end));    
%     s_ms(i,:,:)  = fluid.saturation(states_cm(end));
%     for j = 1:G.cells.num
%         kr_ms(i,j,:) = fluid.relperm(s_ms(i,j,:), states_cm(end));
%         lam_w_ms(i,j)= kr_ms(i,j,1) / mu_ms(i,1); % water mobility
%         lam_o_ms(i,j)= kr_ms(i,j,2) / mu_ms(i,2); % oil mobility
%     end
%     F_w_ms(i,:) = lam_w_ms(i,:) ./ (lam_w_ms(i,:) + lam_o_ms(i,:));
%     PVI_ms(i,1) = W_ms(1).val(1) * (dt*i) / sum(pvCart.*rock.poro);
end
p_f_mean     = zeros(4,1);
p_f_mean(1)  = mean(p_L2);
p_f_mean(2)  = mean(p_Linf);
p_f_mean(3)  = mean(f_L2);
p_f_mean(4)  = mean(f_Linf);
p_f_mean     = p_f_mean';

s_f_mean     = zeros(3,1);
s_f_mean(1)  = mean(s_L1(2:end));
s_f_mean(2)  = mean(s_L2(2:end));
s_f_mean(3)  = mean(s_Linf(2:end));
s_f_mean     = s_f_mean';
% Fw_ms_prod   = F_w_ms(:,prodC);

% err      = zeros(N_step,1);
% wfp      = zeros(N_step,2);
% for i = 1:N_step
%     err(i,1) = sum(abs(states_cm(i).s(:,1) - states(i).s(:,1)))./sum(states(i).s(:,1));
%     wfp(i,:) = [states(i).s(prodC,1), states_cm(i).s(prodC,1)];    
% end
% names = {'Finescale', 'Direct Mimetic'};
% 
% close all; plotToolbar(G, states);
% axis equal tight off
% daspect([1 1 0.2])
% % view(85, 20);
% % plotWell(G, W);
% title(names{1});
% colorbar('horiz')
% 
% figure; plotToolbar(G, states_cm);
% axis equal tight off
% % daspect([1 1 0.2])
% % view(85, 20);
% % plotWell(G, W);
% title(names{2});
% colorbar('horiz')

fprintf('Simulation is done!\n')
