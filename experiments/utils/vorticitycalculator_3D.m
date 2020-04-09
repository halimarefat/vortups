function vor = vorticitycalculator_3D(G, faceFlux)
% Compute fine-scale vorticity using face-based flux field
%
% SYNOPSIS:
%   vor = vorticitycalculator_3D(G, flux)
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
%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

 vel = faceFlux2cellVelocity(G, faceFlux);
 dims =[G.cartDims(1), G.cartDims(2), G.cartDims(3)];
 dim  = G.cartDims(1)*G.cartDims(2)*G.cartDims(3);
 ut  = zeros(dim,1);
 vt  = zeros(dim,1);
 xt  = zeros(dim,1);
 yt  = zeros(dim,1);
 ut(1:G.cells.num) = vel(:,1);
 vt(1:G.cells.num) = vel(:,2);
 xt(1:G.cells.num) = G.cells.centroids(:,1);
 yt(1:G.cells.num) = G.cells.centroids(:,2);
 u   = reshape(ut(:,1), dims);
 v   = reshape(vt(:,1), dims);
 X   = reshape(xt(:,1), dims);
 Y   = reshape(yt(:,1), dims);
%  u   = u';
%  v   = v';
%  X   = X';
%  Y   = Y';
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
