function [ Ac,ind ] = takeunique( ten )
%extract unique columns from a tensor with symmetric columns.
%improve by mex
%
% © Copyright, Oren Levintal, June 13, 2016.

oldvals=ten.vals;
ten.vals=ten.vals(1,:);
[i,j,~]=tfind(ten);

ij=[i,sort(j,2)]; % columns are sorted j1<=j2<=... below i use the function u2c which assumes j1>=j2>=.. so order of columns need to be reversed
[uniqueij,ind]=unique(ij,'rows');

Au=sptensor(uniqueij(:,1),uniqueij(:,end:-1:2),ind(:)',ten.tsize(1),ten.tsize(2:end)); % columns reversed

Ac=u2c(Au);
ind=Ac.vals(:);
Ac.vals=oldvals(:,ind);

end

