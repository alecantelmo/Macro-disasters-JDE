function [index]=getderivs_c(f,v,order,symparams,fname)
%This function is like getderivs, but only unique nonzero derivatives are
%returned. 
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
%     disp(['Differentiating row no. ' num2str(frow) ' of ' fname ' ...'])
    [derivs(frow,:),uncomp(frow,:)]=compderivs_u(f(frow),v(relevant_v(frow,:)),order);
end

% transform uncompression matrices into indices, and create totalloc.
% totalloc(frow,k) stores the total number of nonzero unique k-order derivatives of f(frow).
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
    fun_name=[fname '_d_c' num2str(k)];
    fid = fopen([fun_name '.m'], 'w');
    fprintf(fid,'%s\n', ['function derivs=' fun_name '(vars,params,index)']);
    fprintf(fid,'%s\n', 'n_s=size(vars,1);');
    for i=1:length(v)
        fprintf(fid,'%s\n', [char(v(i)) '=vars(:,' num2str(i) ');']);
    end
    for i=1:length(symparams)
        fprintf(fid,'%s\n', [char(symparams(i)) '=params(' num2str(i) ');']);
    end
    fprintf(fid,'%s\n', ['full_rows=zeros(' num2str(sum(totalloc(:,k))) ',1);']);
    fprintf(fid,'%s\n', ['full_cols=zeros(' num2str(sum(totalloc(:,k))) ',1);']);
    fprintf(fid,'%s\n', ['full_vals=zeros(n_s,' num2str(sum(totalloc(:,k))) ');']);
    fprintf(fid,'%s\n', ['svec=1:n_s;']);
    for frow=1:n_f
        if totalloc(frow,k)>0
            tempderiv=derivs{frow,k};
            fprintf(fid,'%s\n', ['compressed_deriv=zeros(n_s,' num2str(length(tempderiv)) ');']);
            for i=1:length(tempderiv)
                disp_fun_row('compressed_deriv',tempderiv(i),i,fid);
            end
            fprintf(fid,'%s\n', ['uncompressed_deriv=compressed_deriv(:,index.nnz{' num2str(frow) ',' num2str(k) '});']);
            tempstart=sum(totalloc(1:frow-1,k));
            tempend=sum(totalloc(1:frow,k));
            fprintf(fid,'%s\n', ['full_vals(:,' num2str(tempstart+1) ':' num2str(tempend) ')=uncompressed_deriv;']);

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
            tempv=(1:n_v)';
            takev=tempv(relevant_v(frow,:));
            tempcols=takev(tempcols);
            % find location in compressed matrix
            uniqueind=uniquecols2ind(tempcols,n_v);
            % return to index
            index.loc{frow,k}=uniqueind;
            fprintf(fid,'%s\n', ['full_cols(' num2str(tempstart+1) ':' num2str(tempend) ',:)=index.loc{' num2str(frow) ',' num2str(k) '};']);
            fprintf(fid,'%s\n', ['full_rows(' num2str(tempstart+1) ':' num2str(tempend) ')=' num2str(frow) ';']);
        end
    end
    fprintf(fid,'%s\n', ['derivs=sptensor(full_rows,full_cols,full_vals,' num2str(n_f) ',' num2str(nchoosek(n_v+k-1,k)) ');']);
    fclose(fid);
end


