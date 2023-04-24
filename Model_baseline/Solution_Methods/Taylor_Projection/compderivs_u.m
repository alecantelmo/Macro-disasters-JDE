function [derivs,uncomp]=compderivs_u(f,x,order)
% This is like compderivs except that the uncompression matrix returns a
% column vector of the full derivative matrix where non-unique derivatives
% are zeroed.
%
%
% © Copyright, Oren Levintal, June 13, 2016.

n_x=length(x);
derivs=cell(order,1);
uncomp=cell(order,1);

if n_x==1
    tempderiv=f;
    for k=1:order
        tempderiv=jacobian(tempderiv,x);
        derivs{k}=tempderiv;
        uncomp{k}=speye(length(tempderiv));
    end
else
    tempderiv=jacobian(f,x); 
    tempderiv=tempderiv(:);
    
    nnztempderiv=1-logical(tempderiv==0);
    tempind=find(nnztempderiv);
    countdf=sparse(tempind,ones(length(tempind),1),ones(length(tempind),1),numel(tempderiv),1); % counts nonzero derivatives.
    N1=sparse(tempind,1:sum(countdf),ones(1,sum(countdf)),n_x,sum(countdf));
    uncomp{1}=N1; 
    derivs{1}=tempderiv(countdf==1); 

    for k=2:order
        tempderiv=jacobian(derivs{k-1}, x);
        
        nnztempderiv=1-logical(tempderiv==0);
        [i,j]=find(nnztempderiv);

        countdf_short=sparse(i,j,ones(length(i),1),size(tempderiv,1),size(tempderiv,2));
        countdf=uncomp{k-1}*countdf_short; 
        countdf=countdf(:);
        if nnz(countdf)>0
            [U,W]=create_UW(n_x,k,countdf);
            N=sparse(find(countdf),1:sum(countdf),ones(1,sum(countdf)),n_x^k,sum(countdf));
            tempmat=[U'*(N'*kron(speye(n_x),uncomp{k-1}))];
            [i,j]=find(tempmat');
            tempderiv=tempderiv(i);
            [colM,rowM]=find(N);
            [rowW,colW]=find(W);
            tempuncomp=sparse(size(N,1),size(W,1));
            tempeye=speye(size(W,1));
            tempuncomp(colM,:)=(U*U')*tempeye(rowW,:);
        else
            tempderiv=sym(0);
            tempuncomp=sparse(n_x^k,1);
        end
        derivs{k}=tempderiv;
        uncomp{k}=tempuncomp; 
    end
end

end

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
% This code can be used freely for non commercial purposes, provided that
% it is not altered.
%
% (c) Oren Levintal, December 20, 2013

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
