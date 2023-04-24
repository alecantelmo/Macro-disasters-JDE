function [varout1,varout2]=tkron(B2,B1,varargin)
%C=tkron(B2,B1) calculates C=kron(B2,B1).
%
% © Copyright, Oren Levintal, June 13, 2016.

sB1=size(B1.vals,1);
sB2=size(B2.vals,1);

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
        
l1=B1.tsize(1);
l2=B2.tsize(1);
A=fold(spteye(l1*l2),l1,l2);
if nargin>2
    [varout1,varout2]=contraction2(A,B2,B1,varargin{:});
elseif nargin==2
    [varout1,varout2]=contraction2(A,B2,B1);
end

end