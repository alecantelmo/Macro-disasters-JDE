function [X,Xx,Xxx,Xxxx,Xxxxx,ind]=create_X_tensor(N,x_c0,M2,M3,M3x,W2,W3,unique2,unique3,ind,n_ind,maxload,sum)
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(5,1);
end
n_x=x_c0.tsize(1);

Xxx=[]; Xxxx=[]; Xxxxx=[];
n_s=size(x_c0.vals,1);
temp=sptensor(1);
temp.vals=repmat(temp.vals,n_s,1);
X=vconcat(temp,x_c0);

Ix=spteye(n_x,n_s);
Xx=vconcat(sptensor(1,n_x,n_s),Ix);

if N==1
    Xxx=sptensor(1+n_x,n_x^2,n_s);
elseif N>=2
    [temp,ind{1}]=contraction2(M2,x_c0,x_c0,ind{1},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
    twoW2=multscalar(W2,2);
    twoW2.vals=repmat(twoW2.vals,n_s,1);
    [tempx,ind{2}]=contraction2(twoW2,x_c0,Ix,ind{2},n_ind,maxload,sum);
    Xx=vconcat(Xx,unfold(tempx));
    Xxx=vconcat(sptensor(1+n_x,n_x^2,n_s),unfold(twoW2));
end
if N==2
    Xxxx=sptensor(1+n_x+unique2,n_x^3,n_s);
elseif N>=3
    [temp,ind{3}]=contraction3(M3,x_c0,x_c0,x_c0,ind{3},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
    [tempx,ind{4}]=contraction3(M3x,x_c0,x_c0,Ix,ind{4},n_ind,maxload,sum);
    Xx=vconcat(Xx,unfold(tempx));
    Ix6=multscalar(Ix,6);
    [tempxx,ind{5}]=contraction3(W3,x_c0,Ix,Ix6,ind{5},n_ind,maxload,sum);

    Xxx=vconcat(Xxx,unfold(tempxx));
    sixW3=multscalar(W3,6);
    sixW3.vals=repmat(sixW3.vals,n_s,1);
    Xxxx=vconcat(sptensor(1+n_x+unique2,n_x^3,n_s),unfold(sixW3));
end
if N==3
    Xxxxx=sptensor(1+n_x+unique2+unique3,n_x^4,n_s);
end
