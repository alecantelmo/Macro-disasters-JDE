function [ A ] = ptr2d( A,l1,l2 )
%Convert the pointer to 2-D. If the original pointer was obtained from a
%2-D pointer then rows are assumed to be indexed l2,l1. Otherwis, rows
%are assumed to be indexed l1,l2.
%
% © Copyright, Oren Levintal, June 13, 2016.

if size(A.ptr,2)==1
    if isfield(A,'ptr2d')
        if A.ptr2d(1)~=l1 || A.ptr2d(2)~=l2
            error('dimensions of 2D pointer are wrong')
        else
            A=ptr2d_sub(A,l1,l2);
        end
    else % convert to 1,l1,l2 and move l1 to pointer
        A=ptr2d_sub(A,1,l1*l2); 
        A.ptr2d=intarray([1,l1*l2]); % the size of the 2D pointer
        A=ptr2col(A,1);
        A=foldj(A,1,[l1,l2]);
        A=col2ptr(A,1);
        A=ptr1d(A);
        A=rmfield(A,'ptr2d');
        A=col2ptr(A,1);

    end
end

end

function A=ptr2d_sub(A,l1,l2)

tempptr=zeros(l1,l2+1);
tempptr(:,1:l2)=reshape(A.ptr(1:end-1),l2,l1)';
tempptr(:,end)=[tempptr(2:end,1);A.ptr(end)];
A.ptr=intarray(tempptr);
A.tsize=intarray([l1,l2,A.tsize(2:end)]);

end
