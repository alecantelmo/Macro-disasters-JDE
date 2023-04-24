function [ fxxxxx ] = chain5(fv,fvv,fvvv,fvvvv,fvvvvv,vx,vxx,vxxx,vxxxx,vxxxxx,varargin)
%chain5 calculates a 6-dimensional array of fifth derivatives of the composite
%function f(v(x)) with respect to x, using the fifth order multivariate chain
%rule derived in Levintal, Oren, "Fifth Order Perturbation Solution
%to DSGE Models".
%
%   Input arguments: 
%   fv is the Jacobian matrix of f(v) wrt v.
%   fvv is an array of the second derivatives of f wrt v.
%   fvvv is an array of the third derivatives of f wrt v.
%   fvvvv is an array of the fourth derivatives of f wrt v.
%   fvvvvv is an array of the fifth derivatives of f wrt v.
%   vx is the Jacobian matrix of v(x) wrt x.
%   vxx is an array of the second derivatives of v wrt x.
%   vxxx is an array of the third derivatives of v wrt x.
%   vxxxx is an array of the fourth derivatives of v wrt x.
%   vxxxxx is an array of the fifth derivatives of v wrt x.
%
%   All arguments, except fv and vx, can be reshaped in any form by the
%   reshape.m function. fv and vx must be in matrix form, where the rows of
%   fv correspond to the rows of f(v), and the rows of vx correspond to the
%   rows of v(x).
%
% © Copyright, Oren Levintal, June 13, 2016.


[n_f,n_v]=size(fv);
n_x=size(vx,2);

if isempty(varargin)==0
    OMEGA5=varargin{1};
    OMEGA6=varargin{2};
    OMEGA7=varargin{3};
    OMEGA8=varargin{4};
    OMEGA9=varargin{5};
end

% fvvvvv(vx^5)

if isempty(fvvvvv)
    term1=sparse(1,1);
else
    term1=innerkron(n_f,n_v,fvvvvv,vx,vx,vx,vx,vx);
end

% fvvvv(vx^3*vxx)

term2=innerkron(n_f,n_v,fvvvv,vx,vx,vx,vxx);

if isempty(varargin)
    term2=reshape(full(term2),n_f,n_x,n_x,n_x,n_x,n_x);
    term2_OMEGA5=ipermute(term2,[1,5,6,4,3,2])+ipermute(term2,[1,4,6,5,3,2])+ipermute(term2,[1,3,6,5,4,2])+ipermute(term2,[1,2,6,5,4,3])...
        +ipermute(term2,[1,4,5,6,3,2])+ipermute(term2,[1,3,5,6,4,2])...
        +ipermute(term2,[1,2,5,6,4,3])+ipermute(term2,[1,3,4,6,5,2])+ipermute(term2,[1,2,4,6,5,3])+ipermute(term2,[1,2,3,6,5,4]);
    term2_OMEGA5=reshape(term2_OMEGA5,n_f,n_x^5);
else
   term2_OMEGA5=term2*OMEGA5; 
   clear term2 OMEGA5
end

% fvvv(vx^2*vxxx)

term3=innerkron(n_f,n_v,fvvv,vx,vx,vxxx);

if isempty(varargin)
    term3=reshape(full(term3),n_f,n_x,n_x,n_x,n_x,n_x);
    term3_OMEGA6=ipermute(term3,[1,4,5,6,3,2])+ipermute(term3,[1,3,5,6,4,2])+ipermute(term3,[1,2,5,6,4,3])+ipermute(term3,[1,3,4,6,5,2])...
        +ipermute(term3,[1,2,4,6,5,3])+ipermute(term3,[1,2,3,6,5,4])+ipermute(term3,[1,3,4,5,6,2])...
        +ipermute(term3,[1,2,4,5,6,3])+ipermute(term3,[1,2,3,5,6,4])+ipermute(term3,[1,2,3,4,6,5]);
    term3_OMEGA6=reshape(term3_OMEGA6,n_f,n_x^5);
else
    term3_OMEGA6=term3*OMEGA6;
    clear term3 OMEGA6
end

% fvvv(vx*vxx^2)

term4=innerkron(n_f,n_v,fvvv,vx,vxx,vxx);

if isempty(varargin)
    term4=reshape(full(term4),n_f,n_x,n_x,n_x,n_x,n_x);
    term4_OMEGA7=ipermute(term4,[1,4,5,3,6,2])+ipermute(term4,[1,3,5,4,6,2])+ipermute(term4,[1,2,5,4,6,3])...
        +ipermute(term4,[1,3,4,5,6,2])+ipermute(term4,[1,2,4,5,6,3])+ipermute(term4,[1,2,3,5,6,4])...
        +ipermute(term4,[1,4,5,2,6,3])+ipermute(term4,[1,3,5,2,6,4])+ipermute(term4,[1,2,5,3,6,4])...
        +ipermute(term4,[1,3,4,2,6,5])+ipermute(term4,[1,2,4,3,6,5])+ipermute(term4,[1,2,3,4,6,5])...
        +ipermute(term4,[1,3,4,2,5,6])+ipermute(term4,[1,2,4,3,5,6])+ipermute(term4,[1,2,3,4,5,6]);
    term4_OMEGA7=reshape(term4_OMEGA7,n_f,n_x^5);
else
    term4_OMEGA7=term4*OMEGA7;
    clear term4 OMEGA7
end

% fvv(vx*vxxxx)

term5=innerkron(n_f,n_v,fvv,vx,vxxxx);

if isempty(varargin)
    term5=reshape(full(term5),n_f,n_x,n_x,n_x,n_x,n_x);
    term5_OMEGA8=ipermute(term5,[1,3,4,5,6,2])+ipermute(term5,[1,2,4,5,6,3])+ipermute(term5,[1,2,3,5,6,4])...
        +ipermute(term5,[1,2,3,4,6,5])+ipermute(term5,[1,2,3,4,5,6]);
    term5_OMEGA8=reshape(term5_OMEGA8,n_f,n_x^5);
else
    term5_OMEGA8=term5*OMEGA8;
    clear term5 OMEGA8
end

% fvv(vxx*vxxx)

term6=innerkron(n_f,n_v,fvv,vxx,vxxx);

if isempty(varargin)
    term6=reshape(full(term6),n_f,n_x,n_x,n_x,n_x,n_x);
    term6_OMEGA9=ipermute(term6,[1,3,4,5,2,6])+ipermute(term6,[1,2,4,5,3,6])+ipermute(term6,[1,2,3,5,4,6])+ipermute(term6,[1,2,3,4,5,6])...
        +ipermute(term6,[1,3,4,6,2,5])+ipermute(term6,[1,2,4,6,3,5])+ipermute(term6,[1,2,3,6,4,5])...
        +ipermute(term6,[1,2,5,6,3,4])+ipermute(term6,[1,3,5,6,2,4])+ipermute(term6,[1,4,5,6,2,3]);
    term6_OMEGA9=reshape(term6_OMEGA9,n_f,n_x^5);
else
    term6_OMEGA9=term6*OMEGA9;
    clear term6 OMEGA9
end

% fv(vxxxxx)

if isempty(vxxxxx)
    lastterm=sparse(1,1);
else
    lastterm=innerkron(n_f,n_v,fv,vxxxxx);
end

fxxxxx=term1+term2_OMEGA5+term3_OMEGA6+term4_OMEGA7+term5_OMEGA8+term6_OMEGA9+lastterm;

end

