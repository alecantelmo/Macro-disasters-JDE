function [ derivs_uncompressed,ind ] = uncompressderivs( derivs,order,n_v,nnzmat,ind )
%derivs is a sparse tensor that stores unique derivatives of function f w.r.t v.
%order=derivative order
%n_v=length(v);
%nnzmat= a matrix that stores variables that affect f.
%ind=auxiliary index
%
% © Copyright, Oren Levintal, June 13, 2016.

l=derivs.tsize(1);
if isempty(ind)
    original_derivs=derivs;
    if size(derivs.ptr,2)>1
        error('2D ptr not supported')
    end
    
    n_vals=size(derivs.vals,2);
    derivs.vals=1:n_vals;
    sortind=zeros(1,n_vals*factorial(order));
    newcols=zeros(n_vals*factorial(order),order);
    nz=1;
    newptr=zeros(l+1,1);
    newptr(1)=nz;
    subinds='i1';
    subindsvec='i1(:)';
    for j=2:order
        subinds=[subinds ',i' num2str(j)];
        subindsvec=[subindsvec ',i' num2str(j) '(:)'];
    end
    for i=1:l
        deriv_i=takerows(derivs,i); % take row i
        nnzmati=nnzmat(i,:);
        ni=nnz(nnzmati);
        % convert from compress form to unique form without using c2u
        % function because it takes n_v^order memory, and i want only
        % ni^order.
        eval(['[' subinds ']=ndgrid(find(nnzmati));'])
        eval(['nnzind=sub2ind(repmat(n_v,1,order),' subindsvec ');'])
        nnzind=sparse(nnzind,ones(length(nnzind),1),ones(length(nnzind),1),n_v^order,1);
        [loc,nnzind]=nnz_UW(n_v,order,nnzind,'ind');
        smat=sparse(loc,ones(length(loc),1),nnzind,nchoosek(n_v+order-1,order),1);
        eval(['[' subinds ']=ind2sub(repmat(n_v,1,order),full(smat(deriv_i.cols)));'])
        eval(['allind=[' subindsvec '];'])
        deriv_i.cols=intarray(allind); % this is the unique form
        %%%%
        % change column indices
        deriv_i.cols=intarray(nnzmati(deriv_i.cols));
        deriv_i.tsize=intarray([1,repmat(ni,1,order)]);
        % u2c
        deriv_i=u2c(deriv_i);
        
        [~,W]=create_UW(ni,order);
        deriv_i=contraction1(unfold(deriv_i),spmat2sptensor(W));
        %change column indices back to original
        if order==2
            deriv_i=fold(deriv_i,ni,ni);
        elseif order==3
            deriv_i=fold(deriv_i,ni,ni,ni);
        elseif order==4
            deriv_i=fold(deriv_i,ni,ni,ni,ni);
        else
            error('order should be between 2 and 4')
        end
        temp=find(nnzmati);
        deriv_i.cols=temp(deriv_i.cols);

        n_valsi=size(deriv_i.vals,2);
        sortind(nz:nz+n_valsi-1)=deriv_i.vals;
        newcols(nz:nz+n_valsi-1,:)=deriv_i.cols;
        nz=nz+n_valsi;
        newptr(i+1)=nz;
    end

    ind.sortind=sortind(1:nz-1);
    ind.newcols=intarray(newcols(1:nz-1,:));
    ind.ptr=intarray(newptr);
    derivs=original_derivs;
end

derivs_uncompressed=sptensor;
derivs_uncompressed.vals=derivs.vals(:,ind.sortind);
derivs_uncompressed.cols=ind.newcols;
derivs_uncompressed.ptr=ind.ptr;
derivs_uncompressed.tsize=intarray([l,repmat(n_v,1,order)]);


end

