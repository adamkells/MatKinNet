A=zeros(length(unique(data)));
val=1;
x=unique(data);
keyboard
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

%A=sparse(A);

% for i=1:length(A)
%     if A(i,i)==1
%         if sum(A(i,:))==1
%             A(i,:)=[];
%             A(:,i)=[];
%         end
%     end
% end

