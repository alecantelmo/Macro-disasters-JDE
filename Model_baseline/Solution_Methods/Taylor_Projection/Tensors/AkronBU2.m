function [varout1,varout2]=AkronBU2(A,B,varargin)
%C=AkronBU2(A,B) calculates C=A*kron(B,B)U, where U extracts unique columns.
%[C,ind]=AkronBU2(A,B,ind,1,maxload) calculates C using one precomputed
%index ind. If ind=[] the function returns ind. The option maxload controls
%the maximum number of flops performed in parallel by the simd. This is
%relevant only if A or B1,B2 represent sets of s tensors.
%[C,ind]=AkronBU2(A,B,ind,2,maxload) calculates C using two precomputed
%indices ind.
%[C,ind]=AkronBU2(A,B,ind,1,maxload,'sum') sums across s if A
%or B represent sets of s tensors. This option is useful for calculating
%expected value of C.
%[C,ind]=AkronBU2(A,B,ind,2,maxload,'sum') is similar but uses two
%precomputed indices.
%ind=AkronBU2(A,B,1) returns a structure ind with one field
%to speed up calculation.
%ind=AkronBU2(A,B,2) returns a structure ind with two fields. Speedup
%is larger compared to one field, but requires more memory.
%
% © Copyright, Oren Levintal, June 13, 2016.

if nargin<2 || nargin>6
    error('wrong number of input arguments')
end

if A.tsize(2)~=B.tsize(1)
    error('incompatible dimensions')
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
elseif nargin>=4
    maxload=intarray(varargin{3});
end

l=A.tsize(1);
m=A.tsize(2);
if length(A.tsize)~=3
    error('columns must be 2D')
end
if A.tsize(3)~=m
    error('columns are not symmetric')
end
n=B.tsize(2);
if B.colsorted~=1
    error('not sorted')
end
IA=A.ptr;
JA=A.cols;
IB=B.ptr;
JB=B.cols;
if type==0 % calculate A*kron(B,B)U without precomputed index
    varout1=sptensor;
    [IC_1,IC_2] = AkronBU2IC_mex(l,m,n,IA,JA,IB,JB);
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB=B.vals;
    if ~isreal(NA) || ~isreal(NB)
        error('complex numbers not supported')
    end
    [varout1.vals,JC_1] = AkronBU2vecI_mex( l,m,n,IA,JA,NA,IB,...
        JB,NB,IC_1,IC_2,JCrows_1,maxload,dosum);
    varout1.ptr=IC_1;
    varout1.cols=JC_1;
    varout1.tsize=[l,nchoosek(n+1,2)];
elseif type==1 % compute IC
    [varout1.IC_1,varout1.IC_2] = AkronBU2IC_mex(l,m,n,IA,JA,IB,JB);
elseif type==2 % compute IC and JC
    [IC_1,IC_2] = AkronBU2IC_mex(l,m,n,IA,JA,IB,JB);
    JCrows_1=IC_1(end)-1;
    JCrows_2=IC_2(end)-1;
    [ varout1.JC_1,varout1.JC_2 ] = AkronBU2JC_mex( l,m,n,IA,JA,IB,JB,JCrows_1,JCrows_2,IC_2 ); 
    % note: JC_1 has 3 columns. First column stores the location in the
    % compressed matrix. Other columns store the location in a 2-D
    % tensor.
    varout1.IC_1=IC_1;
    varout1.IC_2=IC_2;
elseif type==3 % compute A*kron(B,B)U with one precomputed index IC
    varout1=sptensor;
    if calculate_ind==1
        [IC_1,IC_2] = AkronBU2IC_mex(l,m,n,IA,JA,IB,JB);
    else
        IC_1=ind.IC_1;
        IC_2=ind.IC_2;
    end
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB=B.vals;
    if ~isreal(NA) || ~isreal(NB)
        error('complex numbers not supported')
    end
    [varout1.vals,JC_1] = AkronBU2vecI_mex( l,m,n,IA,JA,NA,IB,...
        JB,NB,IC_1,IC_2,JCrows_1,maxload,dosum);
    varout1.ptr=IC_1;
    varout1.cols=JC_1;
    varout1.tsize=[l,nchoosek(n+1,2)];
    varout2.IC_1=IC_1;
    varout2.IC_2=IC_2;
elseif type==4 % compute A*kron(B,B)U with two precomputed indices IC,JC
    varout1=sptensor;
    if calculate_ind==1
        [IC_1,IC_2] = AkronBU2IC_mex(l,m,n,IA,JA,IB,JB);
        JCrows_1=IC_1(end)-1;
        JCrows_2=IC_2(end)-1;
        [ JC_1,JC_2 ] = AkronBU2JC_mex( l,m,n,IA,JA,IB,JB,JCrows_1,JCrows_2,IC_2 );
    else
        IC_1=ind.IC_1;
        IC_2=ind.IC_2;
        JC_1=ind.JC_1;
        JC_2=ind.JC_2;
    end
    NA=A.vals;
    NB=B.vals;
    if ~isreal(NA) || ~isreal(NB)
        error('complex numbers not supported')
    end
    varout1.vals = AkronBU2vec_mex( l,m,n,IA,JA,NA,IB,...
        JB,NB,IC_1,JC_1,IC_2,JC_2,maxload,dosum);
    varout1.ptr=IC_1;
    varout1.cols=JC_1(:,1); % note that JC_1 stores 2 columns. First column is the compressed format, and other columns are the corresponding 2-D indices
    varout1.tsize=[l,nchoosek(n+1,2)];
    varout2.IC_1=IC_1;
    varout2.IC_2=IC_2;
    varout2.JC_1=JC_1;
    varout2.JC_2=JC_2;
end

if isfield(A,'ptr2d') % if the pointer of A is 2D varout is also 2D
    varout1.ptr2d=A.ptr2d;
end
end