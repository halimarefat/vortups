function p = refineUniformGPU(p, G, indicator, NU, varargin)
%Refine blocks in a partition by uniform partitioning
%
% SYNOPSIS:
%   p = refineUniformGPU(p, G, indicator, NU)
%   p = refineUniformGPU(p, G, indicator, NU, ...
%                            'cartDims',[nx_c ny_c nz_c])
%
% DESCRIPTION:
%   This function refines too large blocks by subdividing the blocks
%   uniformly according to the dimensions given in 'cartDims' (default
%   value is [2,2,2]). Calls partitionUI for the subdivision.
%
% REQUIRED PARAMETERS:
%   p         - Partition vector
%
%   G         - Grid data structure discretising the reservoir model
%               (fine grid, geological model).
%
%   indicator - Cell-wise value of some measure/indicator function used for
%               deciding which blocks to refine.
%
%   NU       - Upper bound
%
% OPTIONAL PARAMETERS:
%
%   cartDims - Dimensions of subdivison of the blocks.
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value dependent upon global verbose
%              settings of function 'mrstVerbose'.
%
% RETURNS:
%   p - Partition vector after refining.

% SEE ALSO:

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
   [p, blockIndicator] = adding(G, p, ijk, blockIndicator, indicator, opt);
end

% Update the partition vector
p = processPartition(G, compressPartition(p));

end
