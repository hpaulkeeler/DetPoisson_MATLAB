% K=funLtoK(L)
% The function funLtoK(L) converts a (non-singular) kernel L matrix into  a (normalized) 
% kernel K matrix. The L matrix has to be semi-positive definite.
%
% Code available here:
% Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_MATLAB

% %TEST
% B=[3, 2, 1; 4, 5,6; 9, 8,7];
% L=B'*B
% L =
%
%   106    98    90
%    98    93    88
%    90    88    86
% K=funLtoK(L)
% K =
%
%     0.7602    0.3348   -0.0906
%     0.3348    0.3320    0.3293
%    -0.0906    0.3293    0.7492

function K=funLtoK(L)
[eigenVectorsLK,eigenValuesL]=eig(L); %eigen decomposition
eigenValuesL=(diag(eigenValuesL)); %eigenvalues of L as vector
eigenValuesK = eigenValuesL./(1+eigenValuesL); %eigenvalues of K
eigenValuesK=diag(eigenValuesK); %%eigenvalues of L as diagonal matrix
K=eigenVectorsLK*eigenValuesK*(eigenVectorsLK'); %recombine from eigen components
K=real(K); %make sure all values are real

end
