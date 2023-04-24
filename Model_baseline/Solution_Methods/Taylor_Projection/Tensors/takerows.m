function [ A,ind ] = takerows( ten,rows )
%extract rows from tensor ten.
%
% © Copyright, Oren Levintal, June 13, 2016.


if nargout==1 %mex version
    A=sptensor;
    IA=ten.ptr;
    JA=ten.cols;
    NA=ten.vals;
    l=ten.tsize(1);
    takel=intarray(rows);
    [IC,JC,NC] = takerows_mex(IA,JA,NA,l,takel);
    newl=length(takel);
    A.ptr=IC(1:newl+1);
    A.cols=JC(1:A.ptr(end)-1,:);
    A.vals=NC(:,1:A.ptr(end)-1);
    A.tsize=ten.tsize;
    A.tsize(1)=newl;
else
    l=length(rows);
    temp=zeros(ten.tsize(1),1);
    temp(rows)=1:l;
    [i,j,vals]=tfind(ten);
    i=temp(i);
    A=sptensor(i(i~=0),j(i~=0,:),vals(:,i~=0),l,ten.tsize(2:end));
    
    ind=find(i);
end

end

