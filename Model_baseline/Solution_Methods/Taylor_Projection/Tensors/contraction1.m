function [varout1,varout2]=contraction1(A,B1,varargin)
%C=contraction1(A,B1) calculates C=A*B1.
%[C,ind]=contraction1(A,B1,ind,1,maxload) calculates C=A*B1 using one precomputed
%index ind. If ind=[] the function returns ind. The option maxload controls
%the maximum number of flops performed in parallel by the simd. This is
%relevant only if A or B1 represent sets of s tensors.
%[C,ind]=contraction1(A,B1,ind,2,maxload) calculates C=A*B1 using two precomputed
%indices ind.
%[C,ind]=contraction1(A,B1,ind,1,maxload,'sum') sums across s if A or B1
%represent sets of s tensors. This option is useful for calculating
%expected value of C.
%[C,ind]=contraction1(A,B1,ind,2,maxload,'sum') is similar but uses two
%precomputed indices.
%ind=contraction1(A,B1,1) returns a structure ind with one field
%to speed up calculation.
%ind=contraction1(A,B1,2) returns a structure ind with two fields. Speedup
%is larger compared to one field, but requires more memory.
%
% © Copyright, Oren Levintal, June 13, 2016.

if nargin<2 || nargin>6
    error('wrong number of input arguments')
end

if A.tsize(2)~=B1.tsize(1)
    error('incompatible dimensions')
end
if size(A.ptr,2)>1
    error('use ptr1d to convert A to 1D')
end
dosum=0;
if nargin==6
    if strcmp(varargin{4},'sum')
        dosum=1;
    elseif ~strcmp(varargin{4},'vec')
        error(['sixth argument must be charater ' '''' 'sum' '''' ' or ' '''' 'vec' ''''])
    end
end
if nargin==2
    type=0; % index is not available
    maxload=intarray(1);
end
if nargin==3
    if varargin{1}==1
        type=1;
    elseif varargin{1}==2
        type=2;
    else
        error('number of precomputed indices should be 1 or 2')
    end
end
calculate_ind=0;
if nargin>=4
    if varargin{2}==1
        type=3;
    elseif varargin{2}==2
        type=4;
    else
        error('number of precomputed indices should be 1 or 2')
    end
    if isempty(varargin{1})
        calculate_ind=1;
    else
        ind=varargin{1};
        if ~isa(ind,'struct')
            error('index variable incorrectly defined')
        end
        if ~isfield(ind,'IC_1')
            error('index variable incorrectly defined')
        end
        if varargin{2}==2
            if ~isfield(ind,'JC_1')
                error('index variable incorrectly defined')
            end
        end
    end
end
if nargin==4
    maxload=intarray(1);
elseif nargin>=5
    maxload=intarray(varargin{3});
end

l=A.tsize(1);
n1=B1.tsize(2);
IA=A.ptr;
JA1=A.cols;
IB1=B1.ptr;
JB1=B1.cols;

ptr2D=0;
if size(B1.ptr,2)>1
    ptr2D=1;
    n1=B1.tsize(3);
end

if type==0 % calculate A*B1 without precomputed index
    varout1=sptensor;
    [ IC_1 ] = contraction1IC_mex( l,n1,IA,JA1,IB1,JB1 );
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB1=B1.vals;
    if ~isreal(NA) || ~isreal(NB1)
        error('complex numbers not supported')
    end
    [varout1.vals,JC1_1] = contraction1vecI_mex( l,n1,IA,JA1,NA,IB1,JB1,NB1,IC_1,JCrows_1,maxload,dosum);
    varout1.ptr=IC_1;
    varout1.cols=JC1_1;
    varout1.tsize=[l,n1];
elseif type==1 % compute IC_1
    if ptr2D==0
        [ varout1.IC_1 ] = contraction1IC_mex( l,n1,IA,JA1,IB1,JB1 );
    else
        [ varout1.IC_1 ] = contraction1IC2D_mex( l,n1,IA,JA1,IB1,JB1 );
    end
elseif type==2 % compute IC_1 and JC_1
    if ptr2D==0
        [ IC_1 ] = contraction1IC_mex( l,n1,IA,JA1,IB1,JB1 );
        JCrows_1=IC_1(end)-1;
        [ JC_1 ] = contraction1JC_mex( l,n1,IA,JA1,IB1,JB1,JCrows_1 );
    else
        [ IC_1 ] = contraction1IC2D_mex( l,n1,IA,JA1,IB1,JB1 );
        JCrows_1=IC_1(end)-1;
        [ JC_1 ] = contraction1JC2D_mex( l,n1,IA,JA1,IB1,JB1,JCrows_1 );
    end
    varout1.IC_1=IC_1;
    varout1.JC_1=JC_1;
elseif type==3 % compute A*B1 with one precomputed index IC
    varout1=sptensor;
    if calculate_ind==1
        if ptr2D==0
            [ IC_1 ] = contraction1IC_mex( l,n1,IA,JA1,IB1,JB1 );
        else
            [ IC_1 ] = contraction1IC2D_mex( l,n1,IA,JA1,IB1,JB1 );
        end
    else
        IC_1=ind.IC_1;
    end
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB1=B1.vals;
    if ~isreal(NA) || ~isreal(NB1)
        error('complex numbers not supported')
    end
    if ptr2D==0
        [varout1.vals,JC_1] = contraction1vecI_mex( l,n1,IA,JA1,NA,IB1,JB1,NB1,IC_1,JCrows_1,maxload,dosum);
        varout1.tsize=[l,n1];
        varout1.ptr=IC_1;
        varout1.cols=JC_1;
        varout2.IC_1=IC_1;
    else
        [varout1.vals,JC_1] = contraction1vecI2D_mex( l,n1,IA,JA1,NA,IB1,JB1,NB1,IC_1,JCrows_1,maxload,dosum);
        varout1.tsize=intarray([l,B1.tsize(2),n1]);
        varout1.ptr=IC_1;
        varout1.cols=JC_1;
        varout1.ptr2d=varout1.tsize(1:2);
        varout2.IC_1=IC_1;
    end
elseif type==4 % compute A*B1 with two precomputed indices IC,JC
    varout1=sptensor;
    if calculate_ind==1
        if ptr2D==0
            [ IC_1 ] = contraction1IC_mex( l,n1,IA,JA1,IB1,JB1 );
            JCrows_1=IC_1(end)-1;
            [ JC_1 ] = contraction1JC_mex( l,n1,IA,JA1,IB1,JB1,JCrows_1 );
        else
            [ IC_1 ] = contraction1IC2D_mex( l,n1,IA,JA1,IB1,JB1 );
            JCrows_1=IC_1(end)-1;
            [ JC_1 ] = contraction1JC2D_mex( l,n1,IA,JA1,IB1,JB1,JCrows_1 );
        end
    else
        IC_1=ind.IC_1;
        JC_1=ind.JC_1;
    end
    NA=A.vals;
    NB1=B1.vals;
    if ~isreal(NA) || ~isreal(NB1)
        error('complex numbers not supported')
    end
    if ptr2D==0
        varout1.vals = contraction1vec_mex( l,n1,IA,JA1,NA,IB1,JB1,NB1,IC_1,JC_1,maxload,dosum);
        varout1.ptr=IC_1;
        varout1.cols=JC_1;
        varout1.tsize=[l,n1];
        varout2.IC_1=IC_1;
        varout2.JC_1=JC_1;
    else
        [varout1.vals] = contraction1vec2D_mex( l,n1,IA,JA1,NA,IB1,JB1,NB1,IC_1,JC_1,maxload,dosum);
        varout1.tsize=intarray([l,B1.tsize(2),n1]);
        varout1.ptr=IC_1;
        varout1.cols=JC_1;
        varout1.ptr2d=varout1.tsize(1:2);
        varout2.IC_1=IC_1;
        varout2.JC_1=JC_1;
    end
end

if isfield(A,'ptr2d') % if the pointer of A is 2D varout is also 2D
    if isfield(varout1,'ptr2d')
        varout1.ptr2d=[A.ptr2d,varout1.ptr2d(2)];
    else
        varout1.ptr2d=[A.ptr2d];
    end
end

end