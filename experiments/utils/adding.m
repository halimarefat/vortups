function [p, blockIndicator] = adding(G, p_in, ijk, blockIndicator_in, indicator, opt)
   [v, block]   = max(blockIndicator_in);
   cellsInBlock = find(p_in==block);
   
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
   part = part + max(p_in);
   p_in(cellsInBlock) = part;
   p = p_in;

   % Updating indicator value for the blocks
   blockIndicator = accumarray(p, indicator);
end 