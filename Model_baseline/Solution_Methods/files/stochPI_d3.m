function derivs=stochPI_d3(vars,params,index)
n_s=size(vars,1);
logcp=vars(:,1);
logc=vars(:,2);
logkp=vars(:,3);
logap=vars(:,4);
logk=vars(:,5);
loga=vars(:,6);
mpkp=vars(:,7);
logmpkp=vars(:,8);
mp=vars(:,9);
logmp=vars(:,10);
logoutput=vars(:,11);
output=vars(:,12);
k=vars(:,13);
c=vars(:,14);
kp=vars(:,15);
invk=vars(:,16);
BETA=params(1);
GAMMA=params(2);
ALPHA=params(3);
RHO=params(4);
DELTA=params(5);
SIGMA=params(6);
full_rows=zeros(2,1);
full_cols=zeros(2,3);
full_vals=zeros(n_s,2);
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logmpkp);
uncompressed_deriv=compressed_deriv(:,index.nnz{7,3});
full_vals(:,1:1)=uncompressed_deriv;
full_cols(1:1,:)=index.loc{7,3};
full_rows(1:1)=7;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logmp);
uncompressed_deriv=compressed_deriv(:,index.nnz{9,3});
full_vals(:,2:2)=uncompressed_deriv;
full_cols(2:2,:)=index.loc{9,3};
full_rows(2:2)=9;
derivs=sptensor(full_rows,full_cols,full_vals,16,repmat(16,1,3));
