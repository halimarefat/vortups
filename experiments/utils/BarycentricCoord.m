clc;
clear;
close all;

T = [10; 10; 100];
P = [0,0; 
     2,0; 
     1,2;
     0,0];

count = 0;
O     = zeros(121,2); 
T_O   = zeros(121,1);
for i = 0:0.2:2
    for j = 0:0.2:2
        count = count + 1;
        O(count, 1) = i;
        O(count, 2) = j;
    end
end

% O = [0.2,1.2];
lambda_1 = zeros(121,1);
lambda_2 = zeros(121,1);
lambda_3 = zeros(121,1);

det_T    = ( (P(1,2) - P(2,2)) * (P(3,1) - P(2,1)) + (P(2,1) - P(1,1)) * (P(3,2) - P(2,2)) );
for i = 1:121
    lambda_1(i) = ( (P(1,2) - P(2,2)) * (O(i,1) - P(2,1)) + (P(2,1) - P(3,1)) * (O(i,2) - P(2,2)) ) / det_T; 
    lambda_2(i) = ( (P(2,2) - P(3,2)) * (O(i,1) - P(2,1)) + (P(3,1) - P(2,1)) * (O(i,2) - P(2,2)) ) / det_T; 
    lambda_3(i) = 1 - lambda_1(i) - lambda_2(i);
    T_O(i)      = (lambda_1(i) * T(1,1) + lambda_2(i) * T(2,1) + lambda_3(i) * T(3,1)) ...
                / (lambda_1(i) + lambda_2(i) + lambda_3(i));
end
T_O_resh = reshape(T_O, [11 11]);
xx = repmat(0:0.2:2,11, 1);
yy = repmat((0:0.2:2)',1, 11);
contourf(xx, yy, T_O_resh); hold on
plot(P(:,1),P(:,2), 'k'); hold on
text(P(3,1),P(3,2), int2str(T(1,1)), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'center');
text(P(1,1),P(1,2), int2str(T(2,1)), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'left');
text(P(2,1),P(2,2), int2str(T(3,1)), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'right');