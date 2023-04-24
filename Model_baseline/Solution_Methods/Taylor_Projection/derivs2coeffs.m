function [ coeffs ] = derivs2coeffs( model,G0,varargin )
%The function converts derivatives into a column of unique Taylor
%coefficients.
%
% © Copyright, Oren Levintal, June 13, 2016.

order=length(varargin);
coeffs=G0(:);
n_y=length(coeffs);

if order>=1
    temp=reshape(varargin{1},n_y,[]);
    n_x=size(temp,2);
    coeffs=[coeffs,temp];
end

for k=2:order
    U=model.U{k};
    temp=reshape(varargin{k},n_y,[])*U/factorial(k);
    coeffs=[coeffs,temp];
end

coeffs=coeffs(:);

end

