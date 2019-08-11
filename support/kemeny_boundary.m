function [kemenyR,R,P_EQ,eigs]=kemeny_boundary(K,INV_K,eq,A,jj,param)

% computing hummer-szabo reduced matrix
if jj==0
    %keyboard
    [R,~,~]=hummer_szabo_clustering_A(K,INV_K, eq, A);
elseif jj==1
    R=localeq(K,eq,A);
end
%keyboard
% analyse eigenvalues vectors of clustered matrix R
[Reigs,P_EQ,~,~,~]=spec_decomp(R);
%keyboard
% compute variational parameter
kemenyR=sum(-1./Reigs(2:param+2));
eigs = (-1./Reigs(2:end));

%keyboard
%kemenyR=sum(-1./Reigs(2:(end-1)));
end