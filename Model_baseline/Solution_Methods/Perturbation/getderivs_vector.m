function [index]=getderivs_vector(f,v,order,symparams,fname)
%This function differentiates f w.r.t v up to order, and generates
%vectorized functions that calculate these derivatives. The derivatives are
%stored in a sparse column vector of length n_f*n_v^order, where n_f
%is the size of f and n_v is the size of v.
%
% © Copyright, Oren Levintal, June 13, 2016.


if order<1
    error('order must be at least 1')
end

f=f(:);
v=v(:);

n_f=length(f);
n_v=length(v);

uncomp=cell(n_f,order);
derivs=cell(n_f,order);

% differentiate rows of f separately w.r.t to relevant variables ONLY
relevant_v=logical(jacobian(f(:),v)~=0);
for frow=1:n_f
    [derivs(frow,:),uncomp(frow,:)]=compderivs(f(frow),v(relevant_v(frow,:)),order);
end

% transform uncompression matrices into indices, and create totalloc.
% totalloc(frow,k) stores the total number of nonzero k-order derivatives of f(frow), including all (nonzero) mixed derivatives.
index=[];
totalloc=zeros(n_f,order);
for k=1:order
    for frow=1:n_f
        [index.loc{frow,k},index.nnz{frow,k}]=find(uncomp{frow,k}); 
        totalloc(frow,k)=length(index.loc{frow,k});
        if isempty(index.loc{frow,k})
            index.loc{frow,k}=0;
        end
    end
end

% create matlab functions to calculate the derivatives
for k=1:order
    fun_name=[fname '_d' num2str(k)];
    fid = fopen([fun_name '.m'], 'w');
    fprintf(fid,'%s\n', ['function derivs=' fun_name '(vars,params,index)']);
    for i=1:length(v)
        fprintf(fid,'%s\n', [char(v(i)) '=vars(' num2str(i) ');']);
    end
    for i=1:length(symparams)
        fprintf(fid,'%s\n', [char(symparams(i)) '=params(' num2str(i) ');']);
    end
    fprintf(fid,'%s\n', ['full_rows=zeros(' num2str(sum(totalloc(:,k))) ',1);']);
    fprintf(fid,'%s\n', ['full_vals=zeros(' num2str(sum(totalloc(:,k))) ',1);']);    

    for frow=1:n_f
        if totalloc(frow,k)>0
            tempderiv=derivs{frow,k};
            fprintf(fid,'%s\n', ['compressed_deriv=zeros(' num2str(length(tempderiv)) ',1);']);    
            for i=1:length(tempderiv)
                disp_fun('compressed_deriv',tempderiv(i),i,fid);
            end
            fprintf(fid,'%s\n', ['uncompressed_deriv=compressed_deriv(index.nnz{' num2str(frow) ',' num2str(k) '});']);    
            tempstart=sum(totalloc(1:frow-1,k));
            tempend=sum(totalloc(1:frow,k));
            fprintf(fid,'%s\n', ['full_vals(' num2str(tempstart+1) ':' num2str(tempend) ')=uncompressed_deriv;']);    

            % transform columns to n-dimensions
            tempcols='tempcol1';
            n_relevant_v=sum(relevant_v(frow,:));
            tempdim=num2str(n_relevant_v);
            for tempk=2:k
                tempcols=[tempcols ',tempcol' num2str(tempk)];
                tempdim=[tempdim ',' num2str(n_relevant_v)];
            end
            eval(['[' tempcols ']=ind2sub([' tempdim '],index.loc{' num2str(frow) ',' num2str(k) '});']);
            tempcols=eval(['[' tempcols ']']);
            % each row of f was differentiated w.r.t relevant variables.
            % now, translate dimensions to the full vector of variables.
            tempv=1:n_v;
            takev=tempv(relevant_v(frow,:));
            tempcols=takev(tempcols);
            % create a linear index
            tempdims='tempcols(:,1)';
            for tempi=2:size(tempcols,2)
                tempdims=[tempdims,',tempcols(:,' num2str(tempi) ')'];
            end
            if k>1
                eval(['tempcols=sub2ind(repmat(' num2str(n_v) ',1,' num2str(k) '),' tempdims ');']);
            end
            index.loc{frow,k}=sub2ind([n_f,n_v^k],repmat(frow,length(tempcols),1),tempcols(:));
            
            fprintf(fid,'%s\n', ['full_rows(' num2str(tempstart+1) ':' num2str(tempend) ')=index.loc{' num2str(frow) ',' num2str(k) '};']);    

        end
    end
    fprintf(fid,'%s\n', ['derivs=sparse(full_rows,ones(length(full_vals),1),full_vals,' num2str(n_f*n_v^k) ',1);']);    
    fclose(fid);
end


