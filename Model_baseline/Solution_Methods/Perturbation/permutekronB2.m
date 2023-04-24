
function [ R ] = permutekronB2(order,perms1,n_x,Matrix,varargin )
%permutekronB: This function calculates product of two or three matrices of
%the form:W3*E(kron(hx,zeta,zeta)) or W3*E(kron(hx,zeta,zeta))*(hx*I*I), and
%permutations of this expression, and sums everything.
%
% © Copyright, Oren Levintal, June 13, 2016.

W=varargin{1};
unique=size(W,1);

if nnz(Matrix)==0
    R=sparse(unique,n_x^order);
else
    % 1. Calculate product of first two terms.

    R=W*Matrix;

    % 2. Prepare an index for the permutations.

    if order>1
        ind=reshape([1:n_x^order],repmat(n_x,1,order));
    end

    % 4. Calculate permutations and sum.

    for i=1:size(perms1,1)
        % a. permute the (expected value of the) stochastic part by permutation
        % indices.
        permute_ind=permute(ind,perms1(i,:));
        permute_Matrix=Matrix(permute_ind,permute_ind);

        tempresult=W*permute_Matrix;

        % b. Sum everything.
        R=R+tempresult;
    end

end

