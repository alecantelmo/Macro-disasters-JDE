function derivs=Phi_d_c3(vars,params,index)
n_s=size(vars,1);
logk=vars(:,1);
loga=vars(:,2);
BETA=params(1);
GAMMA=params(2);
ALPHA=params(3);
RHO=params(4);
DELTA=params(5);
SIGMA=params(6);
full_rows=zeros(1,1);
full_cols=zeros(1,1);
full_vals=zeros(n_s,1);
svec=1:n_s;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=0;
uncompressed_deriv=compressed_deriv(:,index.nnz{1,3});
full_vals(:,1:1)=uncompressed_deriv;
full_cols(1:1,:)=index.loc{1,3};
full_rows(1:1)=1;
derivs=sptensor(full_rows,full_cols,full_vals,1,4);
