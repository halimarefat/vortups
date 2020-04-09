function p = myrefineUniform(p, G, indicator, NU, varargin)

if G.griddim==2,
   opt = struct('cartDims', [2, 2], 'verbose',  mrstVerbose);
else
   opt = struct('cartDims', [2, 2, 2], 'verbose',  mrstVerbose);
end
opt = merge_options(opt, varargin{:});

assert(isfield(G, 'cartDims'));

if isfield(G.cells, 'ijkMap'),
   ijk = G.cells.ijkMap;
else
   [ijk{1:G.griddim}] = ind2sub(G.cartDims, G.cells.indexMap(:));
   ijk = [ijk{:}];
end

indicator      = indicator .* G.cells.volumes;
blockIndicator = accumarray(p, indicator);
upper_bound    = NU*sum(indicator)/G.cells.num;

% While there are too large blocks, then refine
while max(blockIndicator) > upper_bound,
   [v, block]   = max(blockIndicator);
   cellsInBlock = find(p==block);

   fprintf('%5f and %4f  \n', max(blockIndicator), block)
   % Build local ijk-map
   ijk_local = ijk(cellsInBlock,:);
   ijk_local = bsxfun(@minus, ijk_local, min(ijk_local)) + 1;

   % Generate a local grid for the current block, such that we can apply
   % partitionUI on it.
   G_loc = extractSubgrid(G, cellsInBlock);
   G_loc.cells.ijkMap = ijk_local;
   G_loc.cartDims     = max(ijk_local);

   if G_loc.griddim == 2,
      indexMap=sub2ind(G_loc.cartDims, ijk_local(:,1), ijk_local(:,2));
   else
      indexMap=sub2ind(G_loc.cartDims, ...
         ijk_local(:,1), ijk_local(:,2), ijk_local(:,3));
   end

   G_loc.cells.indexMap = indexMap;
   % G_loc = computeGeometry(G_loc);

   part = partitionUI(G_loc, min([opt.cartDims; G_loc.cartDims]));
   part = part + max(p);
   p(cellsInBlock) = part;

   % Updating indicator value for the blocks
   blockIndicator = accumarray(p, indicator);
end

% Update the partition vector
p = processPartition(G, compressPartition(p));

end
