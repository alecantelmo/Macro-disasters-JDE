function [ fxx ] = chain2(fv,fvv,vx,vxx)
%chain2 calculates a 3-dimensional array of second derivatives of the composite
%function f(v(x)) with respect to x, using the second order multivariate chain
%rule derived in Levintal, Oren, "Fifth Order Perturbation Solutions
%to DSGE Models".
%
%   Input arguments: 
%   fv is the Jacobian matrix of f(v) wrt v.
%   fvv is an array of the second derivatives of f wrt v.
%   vx is the Jacobian matrix of v(x) wrt x.
%   vxx is an array of the second derivatives of v wrt x.
%
%   All arguments, except fv and vx, can be reshaped in any form by the
%   reshape.m function. fv and vx must be in matrix form, where the rows of
%   fv correspond to the rows of f(v), and the row of vx correspond to the
%   rows of v(x).
%
% © Copyright, Oren Levintal, June 13, 2016.


[n_f,n_v]=size(fv);
n_x=size(vx,2);

term1=innerkron(n_f,n_v,fvv,vx,vx);
term2=fv*reshape(vxx,n_v,n_x^2);

fxx=term1+term2;


end

