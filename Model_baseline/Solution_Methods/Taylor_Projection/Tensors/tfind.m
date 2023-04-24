function [i,j,vals]=tfind(A)
% similar to the function find for sparse matrices.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(A.vals)
    i=zeros(0,1);
    j=zeros(0,length(A.tsize)-1);
else
    if size(A.ptr,2)==1
        x=zeros(size(A.cols,1)+1,1);
        rows=1:A.tsize(1);
        x(A.ptr(1:end-1))=rows;
        x=x(1:end-1);
        x(x>0)=[x(1);diff(x(x>0))];
        x=cumsum(x);
        i=x;
        j=A.cols;
    else
%convert the 2-D pointer to 1-D
        tempptr=A.ptr(:,1:end-1);
        tempptr=tempptr';
        tempptr=[tempptr(:);A.ptr(end)];
        tempsize=[prod(A.tsize(1:2)),A.tsize(3:end)];
        x=zeros(size(A.cols,1)+1,1);
        rows=1:tempsize(1);
        x(tempptr(1:end-1))=rows;
        x=x(1:end-1);
        x(x>0)=[x(1);diff(x(x>0))];
        x=cumsum(x);
        i=x;
        [i2,i1]=ind2sub(A.tsize(2:-1:1),i);
        i=i1;
        j=[i2,A.cols];
    end
end
vals=A.vals;
end