function [ R ] = permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,varargin )
%permutekron: This function calculates expressions such as
%(fvvv*kron(Vx0,Vx1,Vx1))*E(kron(hx,zeta,zeta)), or 
%(fvvv*kron(Vx0,Vx1,Vx1))*E(kron(I,zeta,zeta))*(kron(hx,I,I)), which is
%equivalent. The function also calculates permutations of these
%expressions. The permutations are in perms1 and perms2. The output is the sum
%of all permutations.
%
% © Copyright, Oren Levintal, June 13, 2016.
   
if isempty(Matrix) % no stochastic matrix
    stochastic_matrix_exists=0;
    Matrix=1;
else
    stochastic_matrix_exists=1;
end

if nnz(Matrix)==0 % stochastic matrix is all zero
    R=sparse(n_f,n_x^order);
else
    if isempty(perms2)
        perms2_exists=0; % no permutation for nonstochastic part
        perms2_trivial=1;
    elseif nnz(perms2-repmat(1:size(perms2,2),size(perms2,1),1))==0
        perms2_exists=0; % permutation for nonstochastic part are degenerate
        perms2_trivial=1;
    else
        perms2_exists=1;
        perms2_order=size(perms2,2);
    end
    
    % 2. The second term is the stochastic matrix

    R=Matrix;
    stochR=R;
    % 3. Prepare an index for the permutations.

    if order>1 && ~isempty(perms1) 
        tempind=[1:n_x^order];
        ind1=reshape(tempind,repmat(n_x,1,order));
        if perms2_exists==1
            ind2=reshape(tempind,[n_x.^derivs(end:-1:2),1]);
        end
        clear tempind

%         sparseM=speye(n_x^order);
    end
    
    % 4. Calculate permutations and sum.

    for i=1:size(perms1,1)
        % a. Permute
        if isequal(perms1(i,:),1:order)
            perms1_trivial=1;
        else
            permute_ind1=permute(ind1,perms1(i,:));
            perms1_trivial=0;
        end
        if perms2_exists==1
            if isequal(perms2(i,:),1:perms2_order)
                perms2_trivial=1;
            else
%                 permute_ind2=permute(ind2,perms2(i,:));
                ipermute_ind2=ipermute(ind2,perms2(i,:));
                perms2_trivial=0;
            end
        end
        if perms1_trivial==0
            if perms2_trivial==0
%                 permute_stochR=sparseM(:,permute_ind2)*sparseM(:,permute_ind1)'*stochR(:,permute_ind1);%P2*P1'*M*P1
                permute_stochR=stochR(:,permute_ind1);

                permute_stochR=permute_stochR(permute_ind1(ipermute_ind2),:);
%                 max(max(abs(permute_stochR-ipermute_stochR)))
%                 isequal(permute_stochR,ipermute_stochR)
            else
%                 permute_stochR=sparseM(:,permute_ind1)'*stochR(:,permute_ind1);
                permute_stochR=stochR(:,permute_ind1);
                permute_stochR=permute_stochR(permute_ind1,:);
%                 isequal(permute_stochR,ipermute_stochR)
            end
        else
            if perms2_trivial==0
%                 permute_stochR=sparseM(:,permute_ind2)*stochR;
                permute_stochR=stochR(ipermute_ind2,:);
            else
                permute_stochR=stochR;
            end
        end    
        % b. Sum
        R=R+permute_stochR;
    end

  
    % 1. Calculate the nonstochastic part: For example:
    % fvvv*kron(Vx0,Vx1,Vx1).
    if stochastic_matrix_exists==0
        fmat=varargin{1};

        temp=reshape(fmat,numel(fmat)/n_v,n_v);
        for j=1:length(derivs)-2
            temp=Treshape(bigprod(temp,sparse(varargin{j+1})),numel(temp)/n_v^2*n_x^derivs(j+1),n_v);
        end

        j=length(derivs)-1;
        nonstoch=[reshape([temp*sparse(varargin{j+1})]',numel(temp)/n_v*n_x^derivs(j+1)/n_f,n_f)]';
        R=nonstoch*R;
    else
        nnzx=logical(sum(logical(R),2)); 

        fmat=varargin{1};
        temp=reshape(fmat,numel(fmat)/n_v,n_v);
        for j=1:length(derivs)-2
            vmat=sparse(varargin{j+1});
            % choose nonzero columns of vmat
            nnzx=reshape(nnzx,n_x^(order-derivs(j+1)),n_x^derivs(j+1));
            tempind=logical(sum(nnzx,1)); % identify relevant columns
            vmat=vmat(:,tempind); % extract relevant columns
            % add the zero columns
            [temprow,tempcol,tempvals]=find(vmat);
            ztempind=find(tempind);
            vmat=sparse(temprow,ztempind(tempcol),tempvals,size(vmat,1),n_x^derivs(j+1));

            temp=Treshape(bigprod(temp,vmat),numel(temp)/n_v^2*n_x^derivs(j+1),n_v);
            nnzx=nnzx';
        end

        j=length(derivs)-1;
        vmat=sparse(varargin{j+1});
        % choose nonzero columns of vmat
        nnzx=reshape(nnzx,n_x^(order-derivs(j+1)),n_x^derivs(j+1));
        tempind=logical(sum(nnzx,1)); % identify relevant columns
        vmat=vmat(:,tempind); % extract relevant columns
        % add the zero columns
        [temprow,tempcol,tempvals]=find(vmat);
        ztempind=find(tempind);
        vmat=sparse(temprow,ztempind(tempcol),tempvals,size(vmat,1),n_x^derivs(j+1));
        nonstoch=[reshape([bigprod(temp,vmat)]',numel(temp)/n_v*n_x^derivs(j+1)/n_f,n_f)]';
        R=nonstoch*R;        
    end
end

