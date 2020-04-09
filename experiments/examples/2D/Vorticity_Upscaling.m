clc;
close all;
clear;

mrstModule add vortups agglom upscaling coarsegrid incomp spe10 book mimetic;
endDim1 = [0
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
18
18
19
19
20
20
21
21
22
22
23
23
24
24
25
25
26
26
27
27
28
28
29
29
30
30
31
31
32
32
33
33
34
34
35
35
36
36
37
37
38
38
39
39
40
40
41
41
42
42
43
43
44
44
45
45
46
46
47
47
48
48
49
49
50
50
51
51
52
52
53
53
54
54
55
55
56
56
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57
57];
endDim2 = [0
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
19
20
20
21
21
22
22
23
23
24
24
25
25
26
26
27
27
28
28
29
29
30
30
31
31
32
32
33
33
34
34
35
35
36
36
37
37
38
38
39
39
40
40
41
41
42
42
43
43
44
44
45
45
46
46
47
47
48
48
49
49
50
50
51
51
52
52
53
53
54
54
55
55
56
56
57
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
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215
216
217];
CLDim1 = [
    %32
%58
20
51
10
57
53
5
17
42
13
24
27
22
56
47
33
11
35
36
30
34
55
48
54
45
18
40
39
52
2
43
26
38
6
4
50
14
7
8
37
9
28
19
12
16
59
3
41
49
15
1
31
25
23
21
44
46
29];
CLDim2 = [
    %167
%139
145
165
172
140
177
192
206
136
215
212
134
162
201
190
137
78
175
132
158
154
125
219
207
57
131
220
144
130
55
171
189
146
34
176
161
76
116
45
216
174
179
79
160
58
120
47
164
44
149
191
210
166
143
128
155
91
135
56
127
98
126
61
196
151
54
41
26
183
202
133
96
21
27
30
40
22
129
141
25
39
122
159
2
156
66
85
204
17
67
148
184
153
178
197
182
24
59
112
29
106
218
163
118
168
99
170
68
186
152
115
73
150
111
7
92
31
102
72
181
8
81
94
50
84
80
10
200
75
208
185
43
194
52
195
205
108
35
119
110
121
49
157
16
209
104
28
48
142
147
19
103
138
105
42
95
180
33
83
11
3
36
87
64
107
74
124
4
12
46
6
86
198
217
32
20
213
62
109
114
15
89
18
123
187
97
63
173
38
214
53
88
37
71
82
70
211
60
9
13
23
93
14
51
199
117
113
90
188
69
100
101
77
203
193
65
169
5];

% CLDim1 = [
% %     18
% % 13
% 20
% 12
% 19
% 17
% 14
% 11
% 16
% 15
% 53
% 22
% 21
% 33
% 57
% 52
% 35
% 55
% 8
% 10
% 24
% 48
% 42
% 9
% 23
% 7
% 54
% 41
% 34
% 36
% 56
% 6
% 49
% 40
% 38
% 50
% 47
% 46
% 58
% 25
% 26
% 37
% 43
% 39
% 51
% 27
% 45
% 44
% 5
% 31
% 4
% 32
% 3
% 28
% 2
% 29
% 30
% 59
% 60];
% 
% CLDim2 = [
% %     165
% % 173
% 172
% 166
% 138
% 174
% 137
% 192
% 193
% 194
% 164
% 195
% 175
% 171
% 215
% 136
% 196
% 163
% 218
% 145
% 197
% 206
% 162
% 191
% 157
% 198
% 154
% 199
% 176
% 217
% 43
% 148
% 134
% 216
% 156
% 170
% 167
% 139
% 200
% 42
% 56
% 131
% 33
% 129
% 120
% 115
% 111
% 96
% 77
% 76
% 41
% 97
% 147
% 169
% 168
% 32
% 143
% 121
% 159
% 190
% 155
% 213
% 31
% 50
% 186
% 219
% 214
% 51
% 146
% 114
% 110
% 212
% 188
% 44
% 187
% 135
% 65
% 64
% 112
% 149
% 94
% 107
% 113
% 79
% 128
% 124
% 183
% 27
% 207
% 118
% 63
% 73
% 30
% 49
% 52
% 55
% 106
% 109
% 95
% 47
% 62
% 123
% 53
% 185
% 40
% 108
% 74
% 28
% 141
% 158
% 122
% 80
% 93
% 153
% 61
% 98
% 150
% 5
% 4
% 92
% 119
% 104
% 208
% 48
% 102
% 86
% 161
% 17
% 66
% 72
% 116
% 6
% 179
% 3
% 142
% 117
% 184
% 34
% 85
% 18
% 105
% 75
% 205
% 82
% 87
% 57
% 13
% 103
% 10
% 130
% 36
% 209
% 81
% 14
% 12
% 126
% 88
% 21
% 26
% 160
% 58
% 2
% 180
% 211
% 37
% 182
% 178
% 91
% 90
% 29
% 133
% 210
% 15
% 24
% 99
% 89
% 101
% 202
% 70
% 83
% 152
% 100
% 203
% 38
% 201
% 19
% 11
% 35
% 68
% 132
% 127
% 69
% 71
% 125
% 39
% 7
% 25
% 204
% 140
% 54
% 144
% 20
% 45
% 22
% 9
% 78
% 67
% 23
% 8
% 16
% 220
% 189
% 60
% 151
% 177
% 181
% 46
% 59
% 84];

Layer= 56;
dims = [60,220];
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
histogram(rock.perm(:,1));
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
% plot(s1, kr1), legend('kr_1', 'kr_2')
W    = [];
inj  = 1;
prod = G.cells.num; 
for i = 1:30
W = addWell(W, G, rock, i,      ...
            'Type', 'bhp' , 'Val', 500*barsa(), ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W = addWell(W, G, rock, (G.cells.num-i+1),      ...
            'Type', 'bhp' , 'Val', 200*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]);     
end        
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
SimNum= 21; 
err_mn= zeros(SimNum,1);
n     = 20;
% while (n < SimNum)
display(n);
% p     = partitionUI(G, [6,22]);
endDimindx= n;
CLDim1tmp = [ 32;  58; CLDim1(1:endDim1(endDimindx))];
CLDim2tmp = [139; 167; CLDim2(1:endDim2(endDimindx))];
% CLDim1tmp = [ 13;  18; CLDim1(1:endDim1(endDimindx))];
% CLDim2tmp = [165; 173; CLDim2(1:endDim2(endDimindx))];
p  = partitionLiner(G, CLDim1tmp, CLDim2tmp); 
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
    outlineCoarseGrid(G, CG.partition, 'Color', 'w', 'lineWidth', 0.3);
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
%     err(ind,:) = [sum(abs(rc.s(:,1)    - rf.s(:,1))), ...
%                   sum(abs(rfFrf.s(:,1) - rf.s(:,1)))]./sum(rf.s(:,1));
%     err(ind,:) = [sum((rc.s(:,1)    - rf.s(:,1)).^2), ...
%                   sum((rfFrf.s(:,1) - rf.s(:,1)).^2)]./G.cells.num;
    err(ind,:) = [(max(rc.s(:,1))    - min(rc.s(:,1)))    / (max(rf.s(:,1)) - min(rf.s(:,1))), ...
                  (max(rfFrf.s(:,1)) - min(rfFrf.s(:,1))) / (max(rf.s(:,1)) - min(rf.s(:,1)))]; 
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