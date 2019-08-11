clear all
close all
potential_type=2;
% each element of the vector is the number of nodes within each cluster
if potential_type==0
    nodes=[20,20,20]; % length of vector defines number of clusters
    [K,Adj]=erdosrenyi_N(nodes,[0.7,0.005]); 
elseif potential_type==1
    [K,Adj,v]=szabo(30);
elseif potential_type==2
    [K,Adj,v]=linear_pot(200);
    K=K';
end
K=K';

N=size(K,1); % N is the total number of nodes
N2=sqrt(N);
% do spectral decomposition of the rate matrix of the system
[Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K);
kemeny = sum(-1./Keigs(2:end)); % kemeny constant of system
slow_rels = -1./Keigs(2:end); % relaxation processes
G=graph(Adj);
keyboard
save(['testsystem_nonoise_deeper_' num2str(potential_type) '.mat'])