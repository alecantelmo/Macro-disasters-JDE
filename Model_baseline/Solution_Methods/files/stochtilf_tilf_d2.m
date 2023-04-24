function derivs=stochtilf_tilf_d2(vars,params,index)
n_s=size(vars,1);
logcp=vars(:,1);
logc=vars(:,2);
logkp=vars(:,3);
logap=vars(:,4);
mpkp=vars(:,5);
logmpkp=vars(:,6);
mp=vars(:,7);
logmp=vars(:,8);
BETA=params(1);
GAMMA=params(2);
ALPHA=params(3);
RHO=params(4);
DELTA=params(5);
SIGMA=params(6);
full_rows=zeros(2,1);
full_cols=zeros(2,2);
full_vals=zeros(n_s,2);
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{1,2});
full_vals(:,1:2)=uncompressed_deriv;
full_cols(1:2,:)=index.loc{1,2};
full_rows(1:2)=1;
derivs=sptensor(full_rows,full_cols,full_vals,1,repmat(8,1,2));
