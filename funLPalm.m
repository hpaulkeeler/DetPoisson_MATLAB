% function LPalm=funLPalm(L,indexPalm) 
% This function calculates the L matrix for a Palm distribution
% conditioned on points existing at locations of the statespace indexed by 
% indexPalm. The method is discussed in the paper by  Blaszczyszyn and 
% Keeler (in Section B of the appendix), which is based on Proposition 1.2
% in the paper by Borodin and Rains[1], but an equivalent result appears 
% in the paper by Shirai and Takahashi[3]; see Theorem 6.5 and 
% Corolloary 6.6 in [3].
%
% INPUTS:
% L = A square L(-matrix-)kernel, which must be (semi-)positive-definite. 
% indexPalm = an index set for the conditioned point, where all the points 
% of the underlying statespace correspond to the rows (or columns) of L.
%
% LPalm= The (reduced) Palm version of the L matrix, which is a square 
% matrix with dimension of size(L,1)-length(indexPalm).
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018.
%
% References:
% [1] Blaszczyszyn and Keeler, "Determinantal thinning of point processes 
% with network learning applications", 2018.
% [2] Borodin and Rains, "Eynard-Mehta theorem, Schur process, and their 
% Pfaffian analogs", 2005
% [3] Shirai and Takahashi, "Random point fields associated with certain 
% Fredholm determinants I -- fermion, poisson and boson point", 2003.
%
% Code available here:
% Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_MATLAB

function LPalm=funLPalm(L,indexPalm)
sizeL=size(L,1);
if max(indexPalm)>sizeL
    error('The index for the Palm points is larger than matrix L.');
end
indexRemain=setdiff(1:sizeL,indexPalm);%index of remaining points/locations
identityTemp=eye(sizeL); %identity matrix
identityTemp(indexPalm,indexPalm)=0;%some ones set to zero
invLBoro=eye(sizeL)/(identityTemp+L); 
identityTemp=eye(length(indexRemain));
%Borodin and Rains construction -- see equation (24) in [1]
LPalm=identityTemp/(invLBoro(indexRemain,indexRemain))-identityTemp; 
end

