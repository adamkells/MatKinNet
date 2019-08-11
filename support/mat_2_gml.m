% %Extracting edges from gml file graph
% clear all
% fileName = 'polbooks';
% inputfile = fopen(['../Data/SFI/' fileName '.gml']);
% AA=[];
% l=0;
% k=1;
% while 1
%     % Get a line from the input file
%     tline = fgetl(inputfile);
%     % Quit if end of file
%     if ~ischar(tline)
%         break
%     end;
%     nums = regexp(tline,'\d+','match');
%     if length(nums)
%         if l==1
%             l=0;
%             AA(k,2)=str2num(nums{1});
%             k=k+1;
%             continue;
%         end
%         AA(k,1)=str2num(nums{1});
%         l=1;
%     else
%         l=0;
%         continue;
%     end
% end
% 
% keyboard % need to do a little work to define data from AA at this pause
fileName='SFI';

A=zeros(length(unique(data)));
val=1;
x=unique(data);

for i=1:length(x)
    100*i/length(x)
    if max(max(ismember(data,x(i))))==1
        data(data==x(i))=val;
        val=val+1;
    end
end

for i=1:size(data,1)
    A(data(i,1),data(i,2))=1;
end
A=A+A';

nodes = max(max(data));
for i = 1:nodes
    for j = 1:nodes
        K(i,j)=A(i,j)/sum(A(:,j));
    end
end
for i = 1:nodes
    K(i,i)=0;
    K(i,i)=-sum(K(i,:));
end

K=K';

Adj = A;

N=size(K,1); % N is the total number of nodes
N2=sqrt(N);
% do spectral decomposition of the rate matrix of the system
[Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K);
kemeny = sum(-1./Keigs(2:end)); % kemeny constant of system
slow_rels = -1./Keigs(2:end); % relaxation processes
G=graph(Adj);

save([fileName '.mat'])