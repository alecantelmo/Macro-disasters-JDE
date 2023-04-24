function ten=vconcat(ten1,ten2)
% vertical concatenation of two tensors of same column size.
%
% © Copyright, Oren Levintal, June 13, 2016.

if ~isequal(ten1.tsize(2:end),ten2.tsize(2:end))
    error('incompatible dimensions')
elseif size(ten1.vals,1)~=size(ten2.vals,1)
    error('incompatible states')
else
    ten=sptensor;
    ten.tsize=ten1.tsize;
    ten.tsize(1)=ten1.tsize(1)+ten2.tsize(1);
    ten.ptr=[ten1.ptr;ten2.ptr(2:end)+ten1.ptr(end)-1];
    ten.cols=[ten1.cols;ten2.cols];
    ten.vals=[ten1.vals,ten2.vals];
end


end

