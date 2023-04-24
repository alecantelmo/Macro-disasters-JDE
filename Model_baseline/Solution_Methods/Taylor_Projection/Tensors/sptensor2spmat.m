function [ A ] = sptensor2spmat( ten )
%convert sparse tensor to sparse matrix
%can be improved by mexing 
%
% © Copyright, Oren Levintal, June 13, 2016.

if isfield(ten,'ptr2d') % 2D pointer
    ten=ptr2d(ten,ten.ptr2d(1),ten.ptr2d(2));
    ten=ptr2col(ten,1);
end

[i,j,vals]=tfind(ten);

colsdim=ten.tsize(2:end);
ncolsdim=length(colsdim);
if ncolsdim==1
    cols=j;
elseif ncolsdim==2
    cols=sub2ind(ten.tsize(2:end),j(:,1),j(:,2));
elseif ncolsdim==3
    cols=sub2ind(ten.tsize(2:end),j(:,1),j(:,2),j(:,3));
elseif ncolsdim==4
    cols=sub2ind(ten.tsize(2:end),j(:,1),j(:,2),j(:,3),j(:,4));
else
    error('tensor rank larger than 5')
end

ten.tsize=double(ten.tsize);
A=sparse(double(i),double(cols),vals',ten.tsize(1),prod(ten.tsize(2:end)));

end

