clc;
close all;
clear;

mrstModule add vortups agglom upscaling coarsegrid incomp spe10 book mimetic;

Layer= 59;
[G,~,rock] = getSPE10setup(Layer);
% rock.perm = zeros(G.cells.num,1);
% rock.poro = zeros(G.cells.num,1);
% rock.perm(:,1) = rockSPE10M2.perm(1:G.cells.num,1);
% rock.poro(:,1) = rockSPE10M2.poro(1:G.cells.num,1);
rock.poro(:,1) = max(rock.poro(:,1),1e-3);
% rock.perm(:,1) = 1;
% rock.poro(:,1) = 0.2;

fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
% s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
% plot(s1, kr1), legend('kr_1', 'kr_2')
W    = [];
inj  = 1;
prod = G.cells.num; 
W = addWell(W, G, rock, inj,      ...
            'Type', 'bhp' , 'Val', 500*barsa(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W = addWell(W, G, rock, prod,      ...
            'Type', 'bhp' , 'Val', 200*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]);             
pv    = poreVolume(G, rock);
hT    = computeTrans(G, rock);
Trans = hT2T(G, hT);
S     = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);
rf    = initState(G, W, 500*barsa, [0, 1]);                  
rf    = incompMimetic(rf, G, S, fluid, 'Wells', W);
vel   = faceFlux2cellVelocity(G, rf.flux);
vor   = vorticitycalculator(G, rf.flux);
iVor  = log10(abs(vor));
iVor  = iVor - min(iVor) + 1;
vel  = sqrt(sum(vel.^2,2));
iVel = log10(vel);
iVel = iVel - min(iVel) + 1;
plotCellData(G, iVel, 'EdgeColor', 'none'); axis tight off

pBase = partitionUI(G, [3,6,1]);
p     = refineUniform(pBase, G, iVel, 700);
outlineCoarseGrid(G, p);
CG    = generateCoarseGrid(G,p);
CG    = coarsenGeometry(CG);
plotGrid(CG);

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
rC    = initState(CG, WC, 500*barsa, [0, 1]);                  
rC    = incompMimetic(rC, CG, CS, fluid, 'Wells', WC);

rockFrf.perm(:,1) = rockC.perm(p(:),1);
rockFrf.poro      = rockC.poro(p(:));
hTFrf             = computeTrans(G, rockFrf);
TransFrf          = hT2T(G, hTFrf);
SFrf  = computeMimeticIP(G, rockFrf, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hTFrf]);
rfFrf = initState(G, W, 500*barsa, [0, 1]);                  
rfFrf = incompMimetic(rfFrf, G, SFrf, fluid, 'Wells', W);

T     = 900.0*year;
DT    = T/20;
T_arr = 0:DT:T;

rc    = deal(rf);
wc    = zeros(numel(T_arr), 2); wc (2:end,:) = NaN;
err   = zeros(numel(T_arr), 1); err(2:end,:) = NaN;
t     = DT;
ind   = 1;
rfs_C = zeros(CG.cells.num, 1);
Frfs_C= zeros(CG.cells.num, 1);
title =@(x) title(x,'FontSize',10,'FontWeight','normal');
while (t <= T)
    rf   = implicitTransport(rf   , G, t, rock   , fluid, 'wells', W );
%     rfFrf= implicitTransport(rfFrf, G, t, rockFrf, fluid, 'wells', W );
%     for ix = 1:CG.cells.num
%        ii = find(p(:) == ix);
%        rfs_C (p(ii))  = sum(rf.s(ii(:,1),1))/numel(ii(:,1));
%        Frfs_C (p(ii)) = sum(rfFrf.s(ii(:,1),1))/numel(ii(:,1));
%        ii = 0;
%     end
    rC   = implicitTransport(rC   ,CG, t, rockC  , fluid, 'wells', WC);
    rc.s = rC.s(CG.partition);
    t  = t + DT;
    clf
    axes('position',[.04 .35 .2 .6])
    plotCellData(G, rf.s(:,1), 'EdgeColor', 'none'), axis equal tight off
    title(sprintf('Fine: %d cells', G.cells.num));
    
    axes('position',[.30 .35 .2 .6])
    plotCellData(G, rc.s(:,1), 'EdgeColor', 'none'), axis equal tight off
    outlineCoarseGrid(G, CG.partition, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
    title(sprintf('Coarse: %d cells', CG.cells.num));
    
%     axes('position',[.56 .35 .2 .6])
%     plotCellData(G, rfFrf.s(:,1), 'EdgeColor', 'none'), axis equal tight off
%     outlineCoarseGrid(G, CG.partition, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
%     title(sprintf('Coarse: %d cells', G.cells.num));
    
    ind = ind + 1;
    axes('position',[.06 .05 .42 .25])
    wc(ind,:) = [rf.s(prod, 1), rc.s(prod)]; %, rfFrf.s(prod, 1)];
    plot(T_arr, wc),
    legend('Fine', 'Coarse', 'Frf');
    axis([T_arr(1) T_arr(end) -.05 1.05]), title('Water cut')
    
    axes('position',[.54 .05 .42 .25])
%     err(ind,:) = [sum(abs(rC.s(:,1)   - rfs_C)), ...
%                   sum(abs(Frfs_C      - rfs_C))]; %./sum(rfs_C);
    err(ind,:) = sum(abs(rc.s(:,1)    - rf.s(:,1)))./sum(rf.s(:,1)); %, ...
%                   sum(abs(rfFrf.s(:,1) - rf.s(:,1)))]./sum(rf.s(:,1));
%     err(ind,:) = [sum((rc.s(:,1)    - rf.s(:,1)).^2), ...
%                   sum((rfFrf.s(:,1) - rf.s(:,1)).^2)]./G.cells.num;
    plot(T_arr, err),
    legend('Coarse', 'Frf');
    set(gca,'XLim',[T_arr(1) T_arr(end)]), title('Error')
    drawnow;
    rf    = incompMimetic(rf,    G, S   , fluid, 'Wells', W );
    rfFrf = incompMimetic(rfFrf, G, SFrf, fluid, 'Wells', W );
    rC    = incompMimetic(rC,   CG, CS  , fluid, 'Wells', WC);
end
err_mn(1,1) = mean(err(:,1)) * 100;
% err_mn(1,2) = mean(err(:,2)) * 100;

% end

display('Simulation is done!');