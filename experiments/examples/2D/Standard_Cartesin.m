clc;
clear;
close all;

%% Cartesian fine Geometry is defined and initial solution is done
[G, ~, rock] = getSPE10setup(56);
rock.poro    = max(rock.poro, 1e-4);
W     = [];
inj   = 1;
prod  = G.cells.num;
W     = addWell(W, G, rock, inj,      ...
            'Type', 'rate' , 'Val', 100*meter^3/day(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W     = addWell(W, G, rock, prod,      ...
            'Type', 'bhp' , 'Val', 100*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]); 
MR    = 0.2; % mobility ratio      
fluid = initSimpleFluid('mu' , [   1,MR*1]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);  
hT    = computeTrans(G, rock);
Trans = hT2T(G, hT);
rS    = initState(G, W, 100*barsa, [0, 1]);
S     = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);
rS    = incompMimetic(rS, G, S, fluid, 'Wells', W);
vor   = vorticitycalculator(G, rS.flux);
pv    = poreVolume(G, rock);

Cdims = [6, 7, 1];
p     = partitionUI(G, Cdims);
CG    = generateCoarseGrid(G, p);
CG    = coarsenGeometry(CG);
pvC   = accumarray(CG.partition  ,pv);
rockC = convertRock2Coarse(G, CG, rock);
rockC.perm = upscalePerm(G, CG, rock, 'T', hT, 'S', S);
hTC   = computeTrans(CG, rockC  );
[~,~,WC] = upscaleTrans(CG, hT, 'match_method', 'lsq_flux', ...
                        'bc_method', 'wells', 'wells', W); 
state0 = initState(CG, WC, 100*barsa, [0, 1]);
SC     = computeMimeticIP(CG, rockC, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [CG.cells.faces(:,1), hTC]);
% rSC   = incompMimetic(rSC, CG, SC, fluid, 'Wells', WC);

totTime= 10*year; %1480657.84431263; %
N_step = 100;
dt     = totTime/N_step;

psolve = @(state) incompMimetic(state, CG, SC, fluid, 'Wells', WC);
tsolve = @(state) implicitTransport(state, CG, dt, rockC, fluid, 'wells', WC);

state  = psolve(state0); 
states = state;
mu  = zeros(N_step, 2);
s   = zeros(N_step, CG.cells.num, 2);
kr  = zeros(N_step, CG.cells.num, 2);
lam_w = zeros(N_step, CG.cells.num);
lam_o = zeros(N_step, CG.cells.num);
F_w = zeros(N_step, CG.cells.num);
PVI = zeros(N_step, 1);
for i = 1:N_step
    fprintf('Step %d of %d: ', i, N_step);
    
    fprintf('Solving pressure... ');
    state = psolve(states(end));
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states = [states; tsolve(state)];
    fprintf('Ok!\n');
    mu(i,:)   = fluid.properties(states(end));    
    s(i,:,:)  = fluid.saturation(states(end));
    for j = 1:CG.cells.num
        kr(i,j,:) = fluid.relperm(s(i,j,:), states(end));
        lam_w(i,j)= kr(i,j,1) / mu(i,1); % water mobility
        lam_o(i,j)= kr(i,j,2) / mu(i,2); % oil mobility
    end
    F_w(i,:) = lam_w(i,:) ./ (lam_w(i,:) + lam_o(i,:));
    PVI(i,1) = W(1).val(1) * (dt*i) / sum(pvC.*rockC.poro);
end
Fw_rf_prod = F_w(:,CG.cells.num);
