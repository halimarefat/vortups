function iVor = sphvormap(G, R, B)
% Single-phase vorticity map generator 
%
% SYNOPSIS:
%   iVor = sphvormap(G, R, B)
%
% PARAMETERS:
%   G - 2D grid with Cartesian topology (G, see 'grid_structure') 
%
%   R - Rock parameters.
%
%   B - Flow boundaries such as wells.
%
% RETURNS:
%   vormap - Single-phase vorticity map for defined fine-scale geometry 
%
% REMARK:
%   The computed single-phase vorticity map can be used as an  
%      indicator for refinement using in 2D geometries.  
%
    
    mrstModule add incomp
    locB  = B;
    for i = 1:numel(locB)
        locB(i).compi = 0; 
        if(locB(i).sign == 1) 
            locB(i).compi = 1; 
        end 
    end
    sF  = initSingleFluid('mu' ,    1*centi*poise     , ...
                          'rho', 1014*kilogram/meter^3);
    st0 = initState(G, locB, 100*barsa, 0);
    hT  = computeTrans(G, R);
    S   = computeMimeticIP(G, R, 'InnerProduct', 'ip_quasirt', ...
                        'FaceTrans', [G.cells.faces(:,1), hT]);
    pS  = incompMimetic(st0, G, S, sF, 'Wells', locB);
    vel = faceFlux2cellVelocity(G, pS.flux);
    dim = [G.cartDims(1), G.cartDims(2)];
    [vorarr, ~] = curl(reshape(G.cells.centroids(:,1), dim)', ...
                       reshape(G.cells.centroids(:,2), dim)', ...
                       reshape(vel(:,1), dim)',               ...
                       reshape(vel(:,2), dim)');
    icount = 0;
    vor = zeros(G.cells.num, 1);
    for ix = 1:G.cartDims(2)
        for iy = 1:G.cartDims(1)
            icount = icount + 1;
            vor(icount, 1) = vorarr(ix, iy);
            if (icount > G.cells.num)
                disp('indicies exceed vorticity array size!')
            end 
        end
    end
    iVor = log10(abs(vor));
    iVor = iVor - min(iVor) + 1;    
end 
