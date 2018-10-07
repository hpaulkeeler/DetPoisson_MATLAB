% function indexConfig=funSimSimpleDPP(eigenVectorsL,eigenValuesL)
% This function simulates a determinantal point process (DPP) provided a
% L matrix. It was used to produce the results in the paper by Blaszczyszyn 
% and Keeler[1].
%
% The algorithm exists in Kulesza and Taskar[2]; see Algorithm 1.
% Also see Algorithm 1 in Lavancier, Moller and Rubak[3].
%
% If you use this code in a publication, please cite the paper by
% Blaszczyszyn and Keeler[1].
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.
% [2] Kulesza and Taskar, "Determinantal point processes for machine 
% learning",Now Publisers, 2012
% [3] Lavancier, Moller and Rubak, "Determinantal point process models and
% statistical inference", Journal of the Royal Statistical Society --
% Series B, 2015.
function indexDPP=funSimSimpleDPP(eigenVectorsL,eigenValuesL)

%START - Sampling/simulating DPP - START
%Retrieve eigenvalues and eigenvectors
eigenVectorsK = (eigenVectorsL); 
eigenValuesL=abs(diag(eigenValuesL)); %eigenvalues of L
eigenValuesK = eigenValuesL./(1+eigenValuesL); %eigenvalues of K

%Bernoulli trials (ie coin flips) to determine number of points
indexEigen = find(rand(length(eigenValuesK),1) <= eigenValuesK);

%number of points in the DPP realization
numbPointsDPP=length(indexEigen);
%retrieve eigenvectors corresponding to successful Bernoulli trials
spaceV = eigenVectorsK(:,indexEigen);

indexDPP=zeros(numbPointsDPP,1); %index for final DPP configuration
for ii=numbPointsDPP:-1:1
    %Compute probabilities for each point i
    Prob_i=sum(spaceV.^2,2); %sum across rows
    Prob_i=Prob_i/sum(Prob_i); %normalize
        
    %Choose a point (from 1 to sizeL) using (prob mass function) Prob_i
    indexDPP(ii)=find(rand(1)<=cumsum(Prob_i),1);
    
    %Choose a vector to remove
    jj=find(spaceV(indexDPP(ii),:),1);
    Vj=spaceV(:,jj); %j-th column of V
    spaceV=spaceV(:,[1:jj-1 jj+1:end]);
    
    %Update matrix V by removing Vj component from the space    
    spaceV=spaceV-Vj.*spaceV(indexDPP(ii),:)/Vj(indexDPP(ii));       
    
    %Orthonormalize using Householder method
    [spaceV,~]=qr(spaceV,0);      
end
indexDPP=sort(indexDPP); %sort index in ascending order
%END - Sampling/simulating DPP - END
end



