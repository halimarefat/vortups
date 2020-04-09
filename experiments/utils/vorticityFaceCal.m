function vor = vorticityFaceCal(G)

dims    = [G.cartDims(1,1), G.cartDims(1,2)];
in_indx = zeros(cells-(2*dims(1,1)+2*(dims(1,2)-2)),1);
ou_indx = zeros(2*dims(1,1)+2*(dims(1,2)-2),1);
in_count = 1;
ou_count = 1;
for i = 1:cells
    cell_indx1 = find(N(:,1)==i);
    cell_indx2 = find(N(:,2)==i);
    if (numel(cell_indx2) == 2 && numel(cell_indx1) == 2)
        in_indx(in_count) = i; 
        in_count = in_count + 1;
    else 
        ou_indx(ou_count) = i;
        ou_count = ou_count + 1;
    end
end
neighbours = zeros(cells*4,2);
for i = 1:numel(in_indx(:,1))
    cell_indx1 = find(N(:,1)==in_indx(i,1));
    cell_indx2 = find(N(:,2)==in_indx(i,1));
    neighbours(G.cells.facePos(in_indx(i,1))+0,1) = N(cell_indx2(1),1);
    neighbours(G.cells.facePos(in_indx(i,1))+1,1) = N(cell_indx1(1),2);
    neighbours(G.cells.facePos(in_indx(i,1))+2,1) = N(cell_indx2(2),1);
    neighbours(G.cells.facePos(in_indx(i,1))+3,1) = N(cell_indx1(2),2);
    
    neighbours(G.cells.facePos(in_indx(i,1))+0,2) = 1;
    neighbours(G.cells.facePos(in_indx(i,1))+1,2) = 2;
    neighbours(G.cells.facePos(in_indx(i,1))+2,2) = 3;
    neighbours(G.cells.facePos(in_indx(i,1))+3,2) = 4;
end

for i = 1:numel(ou_indx(:,1))
    cell_indx1 = find(N(:,1)==ou_indx(i,1));
    cell_indx2 = find(N(:,2)==ou_indx(i,1));
    if(isempty(cell_indx2) ~= true)
        if(numel(cell_indx2) == 1)        
            neighbours(G.cells.facePos(ou_indx(i,1))+0,1) = N(cell_indx2(1),1);
            neighbours(G.cells.facePos(ou_indx(i,1))+2,1) = 0;
        end
        if(numel(cell_indx2) == 2)
            neighbours(G.cells.facePos(ou_indx(i,1))+0,1) = N(cell_indx2(1),1);
            neighbours(G.cells.facePos(ou_indx(i,1))+2,1) = N(cell_indx2(2),1);
        end
    else
        neighbours(G.cells.facePos(ou_indx(i,1))+0,1) = 0;
        neighbours(G.cells.facePos(ou_indx(i,1))+2,1) = 0;
    end
    if(isempty(cell_indx1) ~= true) 
        if(numel(cell_indx1) == 1)
            neighbours(G.cells.facePos(ou_indx(i,1))+1,1) = N(cell_indx1(1),2);    
            neighbours(G.cells.facePos(ou_indx(i,1))+3,1) = 0;
        end
        if(numel(cell_indx1) == 2)
            neighbours(G.cells.facePos(ou_indx(i,1))+1,1) = N(cell_indx1(1),2);    
            neighbours(G.cells.facePos(ou_indx(i,1))+3,1) = N(cell_indx1(2),2);
        end
    else
        neighbours(G.cells.facePos(ou_indx(i,1))+1,1) = 0;    
        neighbours(G.cells.facePos(ou_indx(i,1))+3,1) = 0;
    end     
end

for i = 1:cells 
    for j = 0:3        
        if(neighbours(G.cells.facePos(i)+j,1) ~= 0)
            x_dir = G.cells.centroids(i,1) - ...
                G.cells.centroids(neighbours(G.cells.facePos(i)+j,1),1);
            if(x_dir > 0)
                neighbours(G.cells.facePos(i)+j,2) = 1;  
            elseif(x_dir < 0) 
                neighbours(G.cells.facePos(i)+j,2) = 2;             
            end
        end
    end
    for j = 0:3        
        if(neighbours(G.cells.facePos(i)+j,1) ~= 0)
            y_dir = G.cells.centroids(i,2) - ...
                G.cells.centroids(neighbours(G.cells.facePos(i)+j,1),2);
            if(y_dir > 0)
                neighbours(G.cells.facePos(i)+j,2) = 3;  
            elseif(y_dir < 0) 
                neighbours(G.cells.facePos(i)+j,2) = 4;             
            end
        end
    end
end

neighbour_ord = zeros(cells*4,2);

for i = 1:cells     
    neighbour_ord(G.cells.facePos(i)+0,2) = 1;
    neighbour_ord(G.cells.facePos(i)+1,2) = 2;
    neighbour_ord(G.cells.facePos(i)+2,2) = 3;
    neighbour_ord(G.cells.facePos(i)+3,2) = 4;
end
tmp = zeros(4,2);
for i = 1:4
    tmp(i,2) = i; 
end
for i = 1:cells 
    for k = 1:4
        for j = 0:3
            if(neighbours(G.cells.facePos(i)+j,2) == tmp(k,2))
                tmp(neighbours(G.cells.facePos(i)+j,2),1) = neighbours(G.cells.facePos(i)+j,1);
            end
        end
    end    
    for k = 1:4
        neighbour_ord(G.cells.facePos(i)+k-1,1) = tmp(k,1);
    end
    for k = 1:4
        tmp(k,1) = 0; 
    end
end

figure; plotCellData(G,rock.perm(:,1)  ,'EdgeColor','k'); axis tight off
fluid        = initSingleFluid('mu' ,    1*centi*poise     , ...
                               'rho', 1014*kilogram/meter^3);
W = [];
inj  = 1;
prod = cells; 
W = addWell(W, G, rock, inj,      ...
            'Type', 'rate' , 'Val', 1000*meter^3/day, ...
            'Radius', 0.125*meter, 'Sign', 1, ...
            'name', 'I1', 'InnerProduct', 'ip_tpf');
W = addWell(W, G, rock, prod(1),      ...
            'Type', 'bhp' , 'Val', 0*barsa(), ...
            'Radius', 0.125*meter, 'Sign',-1, ...
            'name', 'P1', 'InnerProduct', 'ip_tpf');       
pv    = poreVolume(G, rock);
hT    = computeTrans(G, rock);
Trans = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ hT, [G.faces.num, 1]);
rf    = initState(G, W, 100*barsa);
rf    = incompTPFA(rf, G, Trans, fluid, 'wells', W, 'use_trans', true);

vel   = faceFlux2cellVelocity(G, rf.flux);
vor_f = zeros(cells*4,2);
for i = 1:cells 
    if(neighbour_ord(G.cells.facePos(i)+0,1) ~= 0)
        vor_f(G.cells.facePos(i)+0,1) = vel(i, 1) - vel(neighbour_ord(G.cells.facePos(i)+0,1) ,1)/...
            abs( G.cells.centroids(i,1) - G.cells.centroids(neighbour_ord(G.cells.facePos(i)+0,1) ,1) );
    end
    if(neighbour_ord(G.cells.facePos(i)+1,1) ~= 0)
        vor_f(G.cells.facePos(i)+1,1) = vel(neighbour_ord(G.cells.facePos(i)+1,1) ,1) - vel(i, 1)/...
            abs( G.cells.centroids(neighbour_ord(G.cells.facePos(i)+1,1) ,1) - G.cells.centroids(i,1));
    end
    if(neighbour_ord(G.cells.facePos(i)+2,1) ~= 0)
        vor_f(G.cells.facePos(i)+2,1) = vel(i, 2) - vel(neighbour_ord(G.cells.facePos(i)+2,1) ,2)/...
            abs( G.cells.centroids(i,2) - G.cells.centroids(neighbour_ord(G.cells.facePos(i)+2,1) ,2) );
    end
    if(neighbour_ord(G.cells.facePos(i)+3,1) ~= 0)
        vor_f(G.cells.facePos(i)+3,1) = vel(neighbour_ord(G.cells.facePos(i)+3,1) ,2) - vel(i, 2)/...
            abs( G.cells.centroids(neighbour_ord(G.cells.facePos(i)+3,1) ,2) - G.cells.centroids(i,2));
    end
    vor_f(G.cells.facePos(i)+0,2) = 1;
    vor_f(G.cells.facePos(i)+1,2) = 2;
    vor_f(G.cells.facePos(i)+2,2) = 3;
    vor_f(G.cells.facePos(i)+3,2) = 4;
end
c_faces_ord = zeros(cells*4,3);
tmp = zeros(4,2);
for i = 1:4
    tmp(i,2) = i; 
end
for i = 1:cells 
    for k = 1:4
        for j = 0:3
            if(G.cells.faces(G.cells.facePos(i)+j,2) == tmp(k,2))
                tmp(G.cells.faces(G.cells.facePos(i)+j,2),1) = G.cells.faces(G.cells.facePos(i)+j,1);
            end
        end
    end
    for k = 1:4
        c_faces_ord(G.cells.facePos(i)+k-1,1) = tmp(k,1);
    end
    for k = 1:4
        tmp(k,1) = 0; 
    end
    c_faces_ord(G.cells.facePos(i)+0,3) = 1;
    c_faces_ord(G.cells.facePos(i)+1,3) = 2;
    c_faces_ord(G.cells.facePos(i)+2,3) = 3;
    c_faces_ord(G.cells.facePos(i)+3,3) = 4;
end
for i = 1:cells*4
    c_faces_ord(i,2) = vor_f(i,1);  
end
w_faces = zeros(G.faces.num,2);
for i = 1:G.faces.num
    aux = find(c_faces_ord(:,1) == i);
    w_faces(i,1) = i;
    w_faces(i,2) = c_faces_ord(aux(1),2);
end
w_faces31 = zeros((dims(1,1)+1)*dims(1,2),3);
w_faces32 = zeros(dims(1,1)*(dims(1,2)+1),3);
for i = 1:numel(w_faces31(:,1))
    w_faces31(i,1) = G.faces.centroids(i,1);
    w_faces31(i,2) = G.faces.centroids(i,2);
    w_faces31(i,3) = w_faces(i,2);
end
for i = 1:numel(w_faces32(:,1))
    w_faces32(i,1) = G.faces.centroids(i+numel(w_faces31(:,1)),1);
    w_faces32(i,2) = G.faces.centroids(i+numel(w_faces31(:,1)),2);
    w_faces32(i,3) = w_faces(i+numel(w_faces31(:,1)),2);
end
w_avg_L = zeros(20,1);
for i = 1:20
    w_avg_L(i) = sum(w_faces31(((i-1)*dims(1,1))+(1:dims(1,1)),3))/dims(1,1);
end
faces = zeros(G.faces.num,1);
for i = 1:G.faces.num
    faces(i) = i;
end


figure; plotCellData(cartGrid([dims(1)+1,dims(2)]),-1*abs(w_faces31(:,3))); axis tight off;
figure; plotCellData(cartGrid([dims(1),dims(2)+1]),-1*abs(w_faces32(:,3))); axis tight off;
figure; patch(w_faces32(:,1),w_faces32(:,2),w_faces32(:,3)); axis tight off;
% faces = boundaryFaceIndices(G, 'top', 1:80, 1:20, []);
figure; plotFaces(G, faces(:), 'r'); axis tight off;
display('Vorticity at cell faces are calculated!');
    
end