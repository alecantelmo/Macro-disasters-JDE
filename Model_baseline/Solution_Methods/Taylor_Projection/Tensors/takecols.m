function [ A,ind ] = takecols( ten,cols )
%extract cols from a tensor class ten with symmetric columns.
%improve by mex
%
% © Copyright, Oren Levintal, June 13, 2016.

m=length(cols);
temp=zeros(ten.tsize(2),1);
temp(cols)=1:m;
[i,j,vals]=tfind(ten);
j=temp(j);
prodj=prod(j,2);
A=sptensor(i(prodj~=0),j(prodj~=0,:),vals(:,prodj~=0),ten.tsize(1),repmat(m,1,length(ten.tsize)-1));

ind=find(prodj);

end

