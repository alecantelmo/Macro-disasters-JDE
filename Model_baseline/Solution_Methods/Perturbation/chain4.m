function [ fxxxx ] = chain4(fv,fvv,fvvv,fvvvv,vx,vxx,vxxx,vxxxx,varargin)
%chain4 calculates a 5-dimensional array of fourth derivatives of the composite
%function f(v(x)) with respect to x, using the fourth order multivariate chain
%rule derived in Levintal, Oren, "Fifth Order Perturbation Solution
%to DSGE Models".
%
%   Input arguments: 
%   fv is the Jacobian matrix of f(v) wrt v.
%   fvv is an array of the second derivatives of f wrt v.
%   fvvv is an array of the third derivatives of f wrt v.
%   fvvvv is an array of the fourth derivatives of f wrt v.
%   vx is the Jacobian matrix of v(x) wrt x.
%   vxx is an array of the second derivatives of v wrt x.
%   vxxx is an array of the third derivatives of v wrt x.
%   vxxxx is an array of the fourth derivatives of v wrt x.
%
%   All arguments, except fv and vx, can be reshaped in any form by the
%   reshape.m function. fv and vx must be in matrix form, where the rows of
%   fv correspond to the rows of f(v), and the row of vx correspond to the
%   rows of v(x).
%
% © Copyright, Oren Levintal, June 13, 2016.



[n_f,n_v]=size(fv);
n_x=size(vx,2);

if isempty(varargin)==0
    OMEGA2=varargin{1};
    OMEGA3=varargin{2};
    OMEGA4=varargin{3};
end

if isempty(fvvvv)
    term1=sparse(1,1);
else
    term1=innerkron(n_f,n_v,fvvvv,vx,vx,vx,vx);
end

% fvvv(vx^2*vxx)

term2=innerkron(n_f,n_v,fvvv,vx,vx,vxx);

if isempty(varargin)
    term2=reshape(full(term2),n_f,n_x,n_x,n_x,n_x);
    term2_OMEGA2=ipermute(term2,[1,4,5,3,2])+ipermute(term2,[1,3,5,4,2])+ipermute(term2,[1,2,5,4,3])...
        +ipermute(term2,[1,3,4,5,2])+ipermute(term2,[1,2,4,5,3])+ipermute(term2,[1,2,3,5,4]);
    term2_OMEGA2=reshape(term2_OMEGA2,n_f,n_x^4);
else
    term2_OMEGA2=term2*OMEGA2;
end

% fvv(vx*vxxx)

term3=innerkron(n_f,n_v,fvv,vx,vxxx);

if isempty(varargin)
    term3=reshape(full(term3),n_f,n_x,n_x,n_x,n_x);
    term3_OMEGA3=ipermute(term3,[1,3,4,5,2])+ipermute(term3,[1,2,4,5,3])+ipermute(term3,[1,2,3,5,4])+ipermute(term3,[1,2,3,4,5]);
    term3_OMEGA3=reshape(term3_OMEGA3,n_f,n_x^4);
else
    term3_OMEGA3=term3*OMEGA3;
end

% fvv(vxx^2)

term4=innerkron(n_f,n_v,fvv,vxx,vxx);

if isempty(varargin)
    term4=reshape(full(term4),n_f,n_x,n_x,n_x,n_x);
    term4_OMEGA4=ipermute(term4,[1,3,4,2,5])+ipermute(term4,[1,2,4,3,5])+ipermute(term4,[1,2,3,4,5]);
    term4_OMEGA4=reshape(term4_OMEGA4,n_f,n_x^4);
else
    term4_OMEGA4=term4*OMEGA4;
end

if isempty(vxxxx)
    term5=sparse(1,1);
else
    term5=innerkron(n_f,n_v,fv,vxxxx);
end

fxxxx=term1+term2_OMEGA2+term3_OMEGA3+term4_OMEGA4+term5;

end

