function derivs=prePI_d1(vars,params,index)
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
full_rows=zeros(13,1);
full_cols=zeros(13,1);
full_vals=zeros(n_s,13);
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{1,1});
full_vals(:,1:1)=uncompressed_deriv;
full_cols(1:1,:)=index.loc{1,1};
full_rows(1:1)=1;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{2,1});
full_vals(:,2:2)=uncompressed_deriv;
full_cols(2:2,:)=index.loc{2,1};
full_rows(2:2)=2;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{3,1});
full_vals(:,3:3)=uncompressed_deriv;
full_cols(3:3,:)=index.loc{3,1};
full_rows(3:3)=3;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{4,1});
full_vals(:,4:4)=uncompressed_deriv;
full_cols(4:4,:)=index.loc{4,1};
full_rows(4:4)=4;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{5,1});
full_vals(:,5:5)=uncompressed_deriv;
full_cols(5:5,:)=index.loc{5,1};
full_rows(5:5)=5;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{6,1});
full_vals(:,6:6)=uncompressed_deriv;
full_cols(6:6,:)=index.loc{6,1};
full_rows(6:6)=6;
compressed_deriv=zeros(n_s,2);
compressed_deriv(:,1)=ALPHA;
compressed_deriv(:,2)=1;
uncompressed_deriv=compressed_deriv(:,index.nnz{11,1});
full_vals(:,7:8)=uncompressed_deriv;
full_cols(7:8,:)=index.loc{11,1};
full_rows(7:8)=11;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logoutput);
uncompressed_deriv=compressed_deriv(:,index.nnz{12,1});
full_vals(:,9:9)=uncompressed_deriv;
full_cols(9:9,:)=index.loc{12,1};
full_rows(9:9)=12;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logk);
uncompressed_deriv=compressed_deriv(:,index.nnz{13,1});
full_vals(:,10:10)=uncompressed_deriv;
full_cols(10:10,:)=index.loc{13,1};
full_rows(10:10)=13;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logc);
uncompressed_deriv=compressed_deriv(:,index.nnz{14,1});
full_vals(:,11:11)=uncompressed_deriv;
full_cols(11:11,:)=index.loc{14,1};
full_rows(11:11)=14;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logkp);
uncompressed_deriv=compressed_deriv(:,index.nnz{15,1});
full_vals(:,12:12)=uncompressed_deriv;
full_cols(12:12,:)=index.loc{15,1};
full_rows(12:12)=15;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-exp(-logk);
uncompressed_deriv=compressed_deriv(:,index.nnz{16,1});
full_vals(:,13:13)=uncompressed_deriv;
full_cols(13:13,:)=index.loc{16,1};
full_rows(13:13)=16;
derivs=sptensor(full_rows,full_cols,full_vals,16,repmat(16,1,1));
