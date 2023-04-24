function obj=sptensor(varargin)
% sptensor constructs a sparse tensor.
% T=sptensor(row,cols,vals,rowdim,colsdim) constructs a sparse tensor T.
% The first index of nnz values is stored in vector row. The other indices
% are stored in matrix cols, where size(cols,1)=length(row) and
% size(cols,2) is the number of the other indices. Total number of indices
% (tensor rank) is size(cols,2)+1. The nnz values are stored in vals, where 
% size(vals,2)=size(cols,1)=length(row). It is possible to store a set of s
% tensors with the same sparse structure but different nnz values, by
% storing the nnz values of these tensors as different rows in vals.
% Namely, size(vals,1)=s. 
%
% © Copyright, Oren Levintal, June 13, 2016.

if nargin==0 % empty tensor
    obj.vals=[];
    obj.ptr=[];
    obj.cols=[];
    obj.tsize=[];
elseif nargin==1 % convert a matrix to a sparse tensor. if matrix is sparse better use spmat2sptensor
    T=full(varargin{1});
    if length(size(T))==2
        vals=T';
        obj.vals=vals(:)';
        [m,n]=size(T);
        obj.ptr=(1:n:m*n+1)';
        obj.cols=repmat((1:n)',m,1);
        obj.tsize=[m,n];
        % convert integers to int32 or int64, depending on the
        archstr = computer('arch');
        if strcmp(archstr(end-1:end),'32')
            sys=32;
        elseif strcmp(archstr(end-1:end),'64')
            sys=64;
        else
            error('cannot determine if MATLAB runs on 32-bit or 64-bit platform.')
        end
        if sys==32
            obj.ptr=int32(obj.ptr);
            obj.cols=int32(obj.cols);
            obj.tsize=int32(obj.tsize);
        elseif sys==64
            obj.ptr=int64(obj.ptr);
            obj.cols=int64(obj.cols);
            obj.tsize=int64(obj.tsize);
        end     
    else
        error('full tensor not supported')
    end
elseif nargin==2 || nargin==3 % create a sparse tensor of zeros
    rowdim=varargin{1};
    colsdim=varargin{2};
    if nargin==3
        n_s=varargin{3};
    else
        n_s=1;
    end
    
    obj.vals=zeros(n_s,0);
    obj.cols=zeros(0,length(colsdim));
    obj.ptr=ones(rowdim+1,1);
    obj.tsize=[rowdim,colsdim(:)'];

    % convert integers to int32 or int64, depending on the
    archstr = computer('arch');
    if strcmp(archstr(end-1:end),'32')
        sys=32;
    elseif strcmp(archstr(end-1:end),'64')
        sys=64;
    else
        error('cannot determine if MATLAB runs on 32-bit or 64-bit platform.')
    end
    if sys==32
        obj.ptr=int32(obj.ptr);
        obj.cols=int32(obj.cols);
        obj.tsize=int32(obj.tsize);
    elseif sys==64
        obj.ptr=int64(obj.ptr);
        obj.cols=int64(obj.cols);
        obj.tsize=int64(obj.tsize);
    end     
elseif nargin==5
    row=varargin{1};
    cols=varargin{2};
    vals=varargin{3};
    rowdim=varargin{4};
    colsdim=varargin{5};
    if size(row,1)~=1 && size(row,2)~=1
        error('first argument must be a vector')
    end
    n_vals=size(vals,2);
    if n_vals~=length(row)
        error('length of row index must equal the number of nnz values')
    end
    if size(cols,1)~=n_vals
        error('length of col indices must equal the number of nnz values')
    end
    if ~isa(vals,'double')
        error('nnz values should be of class double')
    end
    % check indices for nonzero tensors
    if n_vals>0
        % check row index
        if ~isequal(round(row),row)
            error('Index into tensor must be an integer.') 
        end
        if max(row)>rowdim
            error('Index exceeds tensor dimensions')
        end
        if min(row)<1
            error('Index into tensor must be positive.')
        end
        % check cols indices
        if size(cols,2)~=length(colsdim)
            error('incompatible column dimensions')
        end
        if ~isequal(round(cols),cols)
            error('Index into tensor must be an integer.') 
        end
        if max(cols(:))>colsdim
            error('Index exceeds tensor dimensions')
        end
        if min(cols(:))<1
            error('Index into tensor must be positive.')
        end
    end
    % convert integers to int32 or int64, depending on the
    archstr = computer('arch');
    if strcmp(archstr(end-1:end),'32')
        sys=32;
    elseif strcmp(archstr(end-1:end),'64')
        sys=64;
    else
        error('cannot determine if MATLAB runs on 32-bit or 64-bit platform.')
    end
    if sys==32
        row=int32(row);
        cols=int32(cols);
        rowdim=int32(rowdim);
        n_vals=int32(n_vals);
    elseif sys==64
        row=int64(row);
        cols=int64(cols);
        rowdim=int64(rowdim);
        n_vals=int64(n_vals);
    end
    % store nnz values in CRS format
%             [ obj.vals,obj.ptr,obj.cols ] = compress( vals,row,cols,rowdim,n_vals );
    [ obj.vals,obj.ptr,obj.cols ] = compress_mex( vals,row,cols,rowdim,n_vals );

    obj.tsize=[rowdim,colsdim(:)'];
else 
    error('wrong number of arguments')
end

end

