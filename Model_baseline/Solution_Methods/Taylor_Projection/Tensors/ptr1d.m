function [ A ] = ptr1d( A )
%convert the 2-D pointer to 1-D
%
% © Copyright, Oren Levintal, June 13, 2016.

if size(A.ptr,2)>1
    tempptr=A.ptr(:,1:end-1);
    tempptr=tempptr';
    A.ptr=[tempptr(:);A.ptr(end)];
    A.tsize=[prod(A.tsize(1:2)),A.tsize(3:end)];
end

end

