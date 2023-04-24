function [varout1,varout2]=tkron3(B3,B2,B1,varargin)
%C=tkron(B3,B2,B1) calculates C=kron(B3,B2,B1).
%
% © Copyright, Oren Levintal, June 13, 2016.

sB1=size(B1.vals,1);
sB2=size(B2.vals,1);
sB3=size(B3.vals,1);

if sB1==1
    if sB2>1
        B1.vals=repmat(B1.vals,sB2,1);
    end
elseif sB1>1
    if sB2==1
        B2.vals=repmat(B2.vals,sB1,1);
    elseif sB2>1 && sB1~=sB2
        error('incompatible states')
    end
end

sB=max(sB1,sB2);

if sB3==1
    if sB>1
        B3.vals=repmat(B3.vals,sB,1);
    end
elseif sB3>1
    if sB==1
        B2.vals=repmat(B2.vals,sB3,1);
        B1.vals=repmat(B1.vals,sB3,1);
    elseif sB>1 && sB3~=sB
        error('incompatible states')
    end
end
        
l1=B1.tsize(1);
l2=B2.tsize(1);
l3=B3.tsize(1);

A=fold(spteye(l1*l2*l3),l1,l2,l3);
if nargin>3
    [varout1,varout2]=contraction3(A,B3,B2,B1,varargin{:});
elseif nargin==3
    [varout1,varout2]=contraction3(A,B3,B2,B1);
end

end