function vor = vorticitycalculator(G, faceFlux)
% Compute fine-scale vorticity using face-based flux field
%
% SYNOPSIS:
%   vor = vorticitycalculator(G, flux)
%
% PARAMETERS:
%   G        - Fine-scale grid with Cartesian topology (G, see 'grid_structure') 
%
%   faceFlux - Vector of fluxes corresponding to face ordering.
%
% RETURNS:
%   vor      - Fine-scale vorticity used as an indicator for refinement 
%
% REMARK:
%  The computed Fine-scale vorticity can be used as an indicator for 
%   refinement using refineUniformShape function.  
%

 vel = faceFlux2cellVelocity(G, faceFlux);
 dim = [G.cartDims(1), G.cartDims(2)];%, G.cartDims(3)]; 
 u   = reshape(vel(:,1), dim);
 v   = reshape(vel(:,2), dim);
 X   = reshape(G.cells.centroids(:,1), dim);
 Y   = reshape(G.cells.centroids(:,2), dim);
 u   = u';
 v   = v';
 X   = X';
 Y   = Y';
 [vorarr, cav] = curl(X, Y, u, v);
 icount = 0;
 vor    = zeros(G.cells.num, 1);
 for ix = 1:G.cartDims(2)
     for iy = 1:G.cartDims(1)
         icount = icount + 1;
         vor(icount, 1) = vorarr(ix, iy);
         if (icount > G.cells.num)
             disp('indicies exceed vorticity array size!')
         end 
     end
 end
end 
