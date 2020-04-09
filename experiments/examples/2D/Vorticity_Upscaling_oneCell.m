clc;
close all;
clear;

mrstModule add vortups agglom upscaling coarsegrid incomp spe10 book mimetic;
endDim1 = [0
1
1
2
2
3
3
4
4
5
5
6
6
7
7
8
8
9
9
10
10
11
11
12
12
13
13
14
14
15
15
16
16
17
17
18
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77];
endDim2 = [0
0
1
1
2
2
3
3
4
4
5
5
6
6
7
7
8
8
9
9
10
10
11
11
12
12
13
13
14
14
15
15
16
16
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17
17];
CLDim1 = [65
71
64
60
63
62
41
78
58
8
72
67
9
21
73
28
66
5
30
48
59
57
61
4
76
37
46
75
2
22
70
36
7
44
53
6
69
51
3
74
1
33
25
24
39
15
45
52
77
20
49
17
10
26
11
14
56
16
12
79
19
34
40
35
47
32
50
68
13
42
43
23
55
18
38
% 27
% 54
31
29];
CLDim2 = [2
5
13
14
11
19
9
17
3
12
18
16
10
6
8
1
% 7
% 15
4];

Layer= 43;
dims = [80,20];
G    = cartGrid(dims);
G    = computeGeometry(G);

[GSPE,~,rockSPE10M2] = getSPE10setup(Layer);
rock.perm = zeros(G.cells.num,1);
rock.poro = zeros(G.cells.num,1);
rock.perm(:,1) = rockSPE10M2.perm(1:G.cells.num,1);
rock.poro(:,1) = rockSPE10M2.poro(1:G.cells.num,1);
rock.poro(:,1) = max(rock.poro(:,1),1e-3);
% rock.perm(:,1) = 1;
% rock.poro(:,1) = 0.2;
% histogram(rock.poro(:,1));
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
% plot(s1, kr1), legend('kr_1', 'kr_2')
W    = [];
inj  = 1;
prod = G.cells.num; 
W = addWell(W, G, rock, inj,      ...
            'Type', 'rate' , 'Val', 1000*meter^3/day(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W = addWell(W, G, rock, prod,     ...
            'Type', 'bhp' , 'Val', 0*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]);          
pv    = poreVolume(G, rock);
hT    = computeTrans(G, rock);
Trans = hT2T(G, hT);
S     = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', ...
                        [G.cells.faces(:,1), hT]);
rf    = initState(G, W, 1000*barsa, [0, 1]);                  
rf    = incompMimetic(rf, G, S, fluid, 'Wells', W);
vel   = faceFlux2cellVelocity(G, rf.flux);
vor   = vorticitycalculator(G, rf.flux);
SimNum= 21; 
err_mn= zeros(SimNum,1);
n     = 1;
% while (n < SimNum)
display(n);
% p     = partitionUI(G, [6,22]);
endDimindx= n;
% CLDim1tmp = [27; 54; CLDim1(1:endDim1(endDimindx))];
% CLDim2tmp = [ 7; 15; CLDim2(1:endDim2(endDimindx))];
CLDim1tmp = [27, 54];
CLDim2tmp = 19;
p  = partitionLiner(G, CLDim1tmp, CLDim2tmp);
% indc  = zeros(G.cells.num,1);
% indc(:,1) = 1;
% icell = find(pbase(:) == 5);
% indc(icell,1) = 2;
% p     = refineUniform(pbase, G, indc, 350, 'cartDims',[25,6]);
plotCellData(G, vor);
outlineCoarseGrid(G, p, 'Color', 'r', 'lineWidth', 1);
CG    = generateCoarseGrid(G,p);
CG    = coarsenGeometry(CG);
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

T     = 1.0*year;
DT    = 0.10*T;
T_arr = 0:DT:T;

rc    = deal(rf);
wc    = zeros(numel(T_arr), 3); wc (2:end,:) = NaN;
err   = zeros(numel(T_arr), 2); err(2:end,:) = NaN;
t     = DT;
ind   = 1;
rfs_C = zeros(CG.cells.num, 1);
Frfs_C= zeros(CG.cells.num, 1);
title =@(x) title(x,'FontSize',10,'FontWeight','normal');
while (t <= T)
    rf   = implicitTransport(rf   , G, t, rock   , fluid, 'wells', W );
    rfFrf= implicitTransport(rfFrf, G, t, rockFrf, fluid, 'wells', W );
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
    outlineCoarseGrid(G, CG.partition, 'Color', 'k', 'lineWidth', 0.3);
    title(sprintf('Coarse: %d cells', CG.cells.num));
    
    ind = ind + 1;
    axes('position',[.06 .05 .42 .25])
    wc(ind,:) = [rf.s(prod, 1), rc.s(prod), rfFrf.s(prod, 1)];
    plot(T_arr, wc),
    legend('Fine', 'Coarse', 'Frf');
    axis([T_arr(1) T_arr(end) -.05 1.05]), title('Water cut')
    
    axes('position',[.54 .05 .42 .25])
%     err(ind,:) = [sum(abs(rC.s(:,1)   - rfs_C)), ...
%                   sum(abs(Frfs_C      - rfs_C))]; %./sum(rfs_C);
    err(ind,:) = [sum(abs(rc.s(:,1)    - rf.s(:,1))), ...
                  sum(abs(rfFrf.s(:,1) - rf.s(:,1)))]./sum(rf.s(:,1));
%     err(ind,:) = [sum((rc.s(:,1)    - rf.s(:,1)).^2), ...
%                   sum((rfFrf.s(:,1) - rf.s(:,1)).^2)]./G.cells.num;
%     err(ind,:) = [(max(rc.s(:,1))    - min(rc.s(:,1)))    / (max(rf.s(:,1)) - min(rf.s(:,1))), ...
%                   (max(rfFrf.s(:,1)) - min(rfFrf.s(:,1))) / (max(rf.s(:,1)) - min(rf.s(:,1)))]; 
    plot(T_arr, err),
    legend('Coarse', 'Frf');
    set(gca,'XLim',[T_arr(1) T_arr(end)]), title('Error')
    drawnow;
    rf    = incompMimetic(rf,    G, S   , fluid, 'Wells', W );
    rfFrf = incompMimetic(rfFrf, G, SFrf, fluid, 'Wells', W );
    rC    = incompMimetic(rC,   CG, CS  , fluid, 'Wells', WC);
end
err_mn(n,1) = mean(err(:,1)) * 100;
err_mn(n,2) = mean(err(:,2)) * 100;
n = n + 1;
% end

display('Simulation is done!');