function p = partitionLiner(G, CLDim1, CLDim2)
% this function partition the given geometry based on desired lines.
% this just for Cartesian grid.

LinDim2 = zeros(G.cartDims(1,2), 3);
LinDim1 = zeros(G.cartDims(1,1), 3);

LinDim2(1,1) = 1;
LinDim2(1,2) = 1;
LinDim2(1,3) = G.cartDims(1,1);
for i = 2:numel(LinDim2(:,1))
    LinDim2(i,1) = i; 
    LinDim2(i,2) = LinDim2(i-1,2) +  G.cartDims(1,1); 
    LinDim2(i,3) = LinDim2(i-1,3) +  G.cartDims(1,1); 
end

LinDim1(1,1) = 1;
LinDim1(1,2) = 1;
LinDim1(1,3) = G.cells.num - G.cartDims(1,1) + 1;

for i = 2:numel(LinDim1(:,1))
    LinDim1(i,1) = i;
    LinDim1(i,2) = LinDim1(i-1,2) + 1;
    LinDim1(i,3) = LinDim1(i-1,3) + 1; 
end

p1 = zeros(G.cells.num,1);
p2 = zeros(G.cells.num,1);
p  = zeros(G.cells.num,1);

k = 1;
for i = 1:numel(LinDim2(:,1))
    for j = LinDim2(i,2):LinDim2(i,3)
        p2(j) = k;        
    end
    if(find(CLDim2(:) == LinDim2(i,1)))
        k = k + 1;
    end
end

k = 1;
for i = 1:numel(LinDim1(:,1))
    for j = LinDim1(i,2):G.cartDims(1,1):LinDim1(i,3)
        p1(j) = k;        
    end
    if(find(CLDim1(:) == LinDim1(i,1)))
        k = k + 1;
    end
end

cc = numel(CLDim1(:,1))+1;
p(1:G.cartDims(1,1)) = p1(1:G.cartDims(1,1));
for i = 2:numel(LinDim2(:,1))
    for j = LinDim2(i,2):LinDim2(i,3)
        p(j) = p(j-G.cartDims(1,1));        
        if( p2(j) ~= p2(j-G.cartDims(1,1)) )                
            p(j) = p(j-G.cartDims(1,1)) + cc; 
        end
    end
end

end