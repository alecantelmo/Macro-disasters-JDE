function [varout1,varout2]=nnz_UW(n,k,nonzero,varargin)
% [U,W]=nnz_UW(n,k,nonzero) creates compression matrices that compress a
% symmetric row vector A, where A(i)~=0 iff nonzero(i)~=0.
% [loc,ind]=nnz_UW(n,k,nonzero,'ind') returns two indices. loc(i) is the
% location in the compressed matrix of elemnt ind(i) in the uncompressed
% matrix.
%
% © Copyright, Oren Levintal, June 13, 2016.

if ~isempty(varargin)
    if strcmp(varargin{1},'sparse')
        dosparse=1;
    elseif strcmp(varargin{1},'ind')
        dosparse=0;
    else
        error('allowed options are sparse or ind')
    end
else
    dosparse=1;
end

I=find(nonzero);

if isempty(I)
    U=sparse(n^k,nchoosek(n+k-1,k));
    W=sparse(nchoosek(n+k-1,k),n^k);
else
    I=I(:);
    origI=I;
    sizeM='n';
    subind='I1';
    for i=2:k
        sizeM=[sizeM ' n'];
        subind=[subind ',I' num2str(i) ];
    end

    eval(['[' subind ']=ind2sub([' sizeM '],I);']); % I is a linear index of all nonzero elements.
    eval(['M=[' subind '];']);

    sortM=sort(M,2); % sort columns
    [M,I]=sortrows(sortM); % sort rows
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
    
    if k>1
        sortM=sortM(unique,:);
        tempmat=repmat(n+k:-1:n+1,n_unique,1);
        loc=(repmat(prod(n+k-1:-1:n),n_unique,1)-prod(tempmat-repmat(sortM(:,1),1,k),2))/factorial(k)+1;

        for j=2:k-1
            tempmat=repmat(n+k-j+1:-1:n+1,n_unique,1);
            loc=loc+(prod(tempmat-repmat(sortM(:,j-1),1,k-j+1),2)-prod(tempmat-repmat(sortM(:,j),1,k-j+1),2))/factorial(k-j+1);
        end
        loc=loc+sortM(:,k)-sortM(:,k-1);
        if dosparse==1
            U= sparse(origI(unique),loc,ones(n_unique,1),n^k,nchoosek(n+k-1,k));
            origI=origI(I); % sort the original index
            W= sparse(loc(group),origI,ones(length(origI),1),nchoosek(n+k-1,k),n^k);
            varout1=U;
            varout2=W;
        else
            varout1=loc;
            varout2=origI(unique);
        end
    elseif k==1
        if dosparse==1
            U=sparse(origI,origI,ones(length(origI),1),n,n);
            W=U;
            varout1=U;
            varout2=W;
        else
            varout1=origI;
            varout2=origI;
        end
    end
end

end


