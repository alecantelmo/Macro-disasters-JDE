function [U,W]=create_UW(n,k,varargin)
% [U,W]=create_UW(n,k) creates two sparse matrices, U and W, that compress and 
% uncompress a symmetric array A with k dimensions and n^k elements. 
% A(:)'*U is a row vector that contains the unique elements of A, 
% and A(:)'=(A(:)'*U)*W. 
% [U,W]=create_UW(n,k,N) creates matrices U and W of a sparse symmetric
% array A with k dimensions and n^k elements. N is a vector of size n^k,
% where N(j)=1 if A(j)~=0 and zero otherwise. If B is a row
% vector with the nonzero elements of the row vector A(:)', then B*U
% is a row vector of the unique nonzero elements of A, and B=(B*U)*W. 
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(varargin)
%     nonzero=ones(n^k,1);
    I=[1:n^k];
else
    nonzero=varargin{1};
    I=find(nonzero);
end

if isempty(I)
    error('The symmetric array is all zero')
end

I=I(:);
% nonzero(nonzero~=0)=1; nonzero=nonzero(:);
sizeM='n';
subind='I1';
for i=2:k
    sizeM=[sizeM ' n'];
    subind=[subind ',I' num2str(i) ];
end

eval(['[' subind ']=ind2sub([' sizeM '],I);']); % I is a linear index of all nonzero elements.
eval(['M=[' subind '];']);

M=sort(M,2); % sort columns
[M,I]=sortrows(M); % sort rows
if size(M,1)==1
    DM=1;
else
    DM=[1;(sum(abs(diff(M)),2)~=0)];
end
group=cumsum(DM);
minI=accumarray(group,I,[],@min);
newi=minI(group);

unique=newi(DM>0); % index of unique nonzero elements

n_unique=length(unique);

U= sparse(unique,[1:n_unique]',ones(n_unique,1),length(I),n_unique);

W= sparse(group,I,ones(length(I),1),n_unique,length(I));

end


