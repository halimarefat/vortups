function node_P = nodePressureCal(node_i, G, state_press, nnc, nav)
%   nnc: nodes_nei_cells
%   nav: nodes_avger

    x_o = G.nodes.coords(node_i,1);
    y_o = G.nodes.coords(node_i,2);
    if (nav(node_i) == 3)
        cells = nnc(node_i,:);
        x     = G.cells.centroids(cells,1);
        y     = G.cells.centroids(cells,2);
        det_T = (y(2) - y(3)) * (x(1) - x(3)) + (x(3) - x(2)) * (y(1) - y(3));
        lam(1)= (y(2) - y(3)) * (x_o  - x(3)) + (x(3) - x(2)) * (y_o  - y(3));
        lam(2)= (y(3) - y(1)) * (x_o  - x(3)) + (x(1) - x(3)) * (y_o  - y(3));
        lam   = lam ./ det_T;
        lam(3)= 1 - lam(1) - lam(2);
        P     = state_press(cells);
        node_P= 0;
        for i = 1:nav(node_i)
            node_P= node_P + lam(i) * P(i);
        end
    else
        cells = nnc(node_i,:);
        indx  = find(cells ~= 0);
        P     = state_press(cells(indx(:)));        
        node_P= sum(P(:))/nav(node_i);    
    end
    
end