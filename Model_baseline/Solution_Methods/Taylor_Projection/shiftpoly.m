function [ G0,varargout ] = shiftpoly( x,c,index,varargin )
%[newG0,newG1,...]=shiftpoly(x,c,index,G0,G1,...) shifts the center of a given 
%polynomial from c to x. The coefficient matrices of the given Polynomial 
%are G0,G1,G2,..... Note that G2,G3,... must have symmetric columns.
%Namely, these are the coefficients of the Taylor series about c.
%The output is the coefficients of the new Polynomial. 
%If index is empty, the variables of the output Polynomial are the same as
%the original Polynomial. To get a Polynomial in a subset of the original
%variables (treating all other vars as constant), specify the index of the
%new variables in index.
%
% © Copyright, Oren Levintal, June 13, 2016.

N=length(varargin)-1;
n_x=length(x);
if isempty(index)
    index=1:n_x;
end

coeffs=cell2mat(varargin);

% create cell
cellX=cell(1,N+1);
m=0;
for i=1:N+1
    newm=m+n_x^(i-1);
    cellX{i}=sparse(size(coeffs,2),n_x^(i-1));
    cellX{i}(m+1:newm,:)=factorial(i-1)*speye(n_x^(i-1));
    m=newm;
end

x_c=x-c;
lastkron=1;
m=1;
for i=1:N % i is the Kronecker power of x-c. All N powers need to be calculated
    newm=m+n_x^i;
    lastkron=kron(lastkron,x_c);
    mj=m;
    for j=1:N+1-i % j-1 is the derivative order. The i Kronecker power appears in derivatives up to order N+1-i
        
        newmj=mj+n_x^(i+j-1);
        if j==1
            A=1;
        else
            A=prod(i+1:i+j-1);
        end
        cellX{j}(mj+1:newmj,:)=A*kron(lastkron,speye(n_x^(j-1)));
        mj=newmj;
    end
    m=newm;
end

xtilx=speye(n_x);
xtilx=xtilx(:,index);

G0=coeffs*cellX{1};

varargout=cell(1,N);
lastkron=1;
for i=1:N
    lastkron=kron(lastkron,xtilx);
    varargout{i}=coeffs*cellX{i+1}*lastkron/factorial(i);

end

end

