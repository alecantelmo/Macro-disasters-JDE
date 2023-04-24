function [varout1,varout2]=AkronBU4i(A,B,ind,n_ind,maxload,dosum,convertind,maxind)
%calculate C=A*kron(B,B,B,B)*U, row by row. Like AkronBU4 but more memory
%efficient.
%[C,ind]=AkronBU4i(A,B,ind,n_ind,maxload,dosum,convertind,maxind).
%convertind is a l-by-m matrix. Each row records the variables that have a
%nonzero effect on A. maxind is the maximum number of such variables across
%rows.
%
% © Copyright, Oren Levintal, June 13, 2016.

if A.tsize(2)~=B.tsize(1)
    error('incompatible dimensions')
end
if strcmp(dosum,'sum')
    dosum=1;
elseif strcmp(dosum,'vec')
    dosum=0;
else
    error(['sixth argument must be charater ' '''' 'sum' '''' ' or ' '''' 'vec' ''''])
end

maxload=intarray(min(1,maxload));
maxind=intarray(maxind);
convertind=intarray(convertind);

if n_ind==1
    type=3;
elseif n_ind==2
    type=4;
else
    error('number of precomputed indices should be 1 or 2')
end

calculate_ind=0;

if isempty(ind)
    calculate_ind=1;
else
    if ~isa(ind,'struct')
        error('index variable incorrectly defined')
    end
    if ~isfield(ind,'IC_1')
        error('index variable incorrectly defined')
    end
    if n_ind==2
        if ~isfield(ind,'JC_1')
            error('index variable incorrectly defined')
        end
    end
end

l=A.tsize(1);
m=A.tsize(2);
if length(A.tsize)~=5
    error('columns must be 4D')
end
if A.tsize(3)~=m || A.tsize(4)~=m || A.tsize(5)~=m
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

if type==3 % compute A*kron(B,B,B,B)*U with one precomputed index IC
    varout1=sptensor;
    if calculate_ind==1
        [IC_1,IC_2,IC_3,IC_4] = AkronBU4iIC_mex(l,m,n,IA,JA,IB,JB,convertind,maxind);
    else
        IC_1=ind.IC_1;
        IC_2=ind.IC_2;
        IC_3=ind.IC_3;
        IC_4=ind.IC_4;
    end
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB=B.vals;
    if ~isreal(NA) || ~isreal(NB)
        error('complex numbers not supported')
    end
    [varout1.vals,JC_1] = AkronBU4ivecI_mex(l,m,n,IA,JA,NA,...
         IB,JB,NB,...
         IC_1,IC_2,IC_3,IC_4,JCrows_1,...
         maxload,dosum,convertind,maxind); 
    varout1.ptr=IC_1;
    varout1.cols=JC_1;
    varout1.tsize=[l,nchoosek(n+3,4)];
    varout2.IC_1=IC_1;
    varout2.IC_2=IC_2;
    varout2.IC_3=IC_3;
    varout2.IC_4=IC_4;
elseif type==4 % compute A*kron(B4,B3,B2,B1) with two precomputed indices IC,JC
    varout1=sptensor;
    if calculate_ind==1
        [IC_1,IC_2,IC_3,IC_4] = AkronBU4iIC_mex(l,m,n,IA,JA,IB,JB,convertind,maxind);
        JCrows_1=IC_1(end)-1;
        JCrows_2=IC_2(end)-1;
        JCrows_3=IC_3(end)-1;
        JCrows_4=IC_4(end)-1;
        [ JC_1,JC_2,JC_3,JC_4 ] = AkronBU4iJC_mex( l,m,n,IA,JA,IB,JB,JCrows_1,JCrows_2,JCrows_3,JCrows_4,IC_2,IC_3,IC_4,convertind,maxind);
    else
        IC_1=ind.IC_1;
        IC_2=ind.IC_2;
        IC_3=ind.IC_3;
        IC_4=ind.IC_4;
        JC_1=ind.JC_1;
        JC_2=ind.JC_2;
        JC_3=ind.JC_3;
        JC_4=ind.JC_4;
    end
    NA=A.vals;
    NB=B.vals;
    if ~isreal(NA) || ~isreal(NB)
        error('complex numbers not supported')
    end
    varout1.vals = AkronBU4ivec_mex(l,m,n,IA,JA,NA,...
         IB,JB,NB,...
         IC_1,JC_1,IC_2,JC_2,IC_3,JC_3,IC_4,JC_4,...
         maxload,dosum,convertind,maxind); 
    varout1.ptr=IC_1;
    varout1.cols=JC_1(:,1);
    varout1.tsize=[l,nchoosek(n+3,4)];
    varout2.IC_1=IC_1;
    varout2.IC_2=IC_2;
    varout2.IC_3=IC_3;
    varout2.IC_4=IC_4;
    varout2.JC_1=JC_1;
    varout2.JC_2=JC_2;
    varout2.JC_3=JC_3;
    varout2.JC_4=JC_4;
end
end