function derivs=stochtilf_tilf_d2(vars,params,index)
n_s=size(vars,1);
auxvar1p=vars(:,1);
nlp=vars(:,2);
logtilqep=vars(:,3);
auxvar1=vars(:,4);
nl=vars(:,5);
logtilqe=vars(:,6);
logqf=vars(:,7);
logtilkstarbackp=vars(:,8);
dp=vars(:,9);
logthetap=vars(:,10);
zap=vars(:,11);
logtilkstarback=vars(:,12);
d=vars(:,13);
logtheta=vars(:,14);
za=vars(:,15);
loghata=vars(:,16);
loghatap=vars(:,17);
loghatz=vars(:,18);
loghatzp=vars(:,19);
u1p=vars(:,20);
u2p=vars(:,21);
logtilvp_tilvss=vars(:,22);
l=vars(:,23);
lp=vars(:,24);
logtilu=vars(:,25);
logtilup=vars(:,26);
loguc=vars(:,27);
logucp=vars(:,28);
logtillambda=vars(:,29);
logtillambdap=vars(:,30);
logtilc=vars(:,31);
logtilcp=vars(:,32);
logmnomp=vars(:,33);
logmden=vars(:,34);
theta=vars(:,35);
thetap=vars(:,36);
tilrp=vars(:,37);
tilq=vars(:,38);
tilqp=vars(:,39);
logtilrp=vars(:,40);
logtilw=vars(:,41);
logtilwp=vars(:,42);
tilcp=vars(:,43);
tilxp=vars(:,44);
u01p=vars(:,45);
u02=vars(:,46);
u03=vars(:,47);
u04=vars(:,48);
logMp=vars(:,49);
u17=vars(:,50);
tilyp=vars(:,51);
wp=vars(:,52);
tildivp=vars(:,53);
MUD=params(1);
RHOTHETA=params(2);
THETABAR=params(3);
SDV_THETA=params(4);
SDV_ZA=params(5);
ALPHA=params(6);
LAMBDA_A=params(7);
GAMMA=params(8);
NU=params(9);
PSI=params(10);
BETA=params(11);
DELTA=params(12);
PHI=params(13);
SCALEPARAM=params(14);
LOGTILUSS=params(15);
full_rows=zeros(31,1);
full_cols=zeros(31,2);
full_vals=zeros(n_s,31);
compressed_deriv=zeros(n_s,2);
compressed_deriv(:,1)=u02.*exp(-u01p.*(GAMMA - 1)).*(GAMMA - 1).^2;
compressed_deriv(:,2)=-exp(-u01p.*(GAMMA - 1)).*(GAMMA - 1);
uncompressed_deriv=compressed_deriv(:,index.nnz{1,2});
full_vals(:,1:3)=uncompressed_deriv;
full_cols(1:3,:)=index.loc{1,2};
full_rows(1:3)=1;
compressed_deriv=zeros(n_s,6);
compressed_deriv(:,1)=u04.*exp(u03);
compressed_deriv(:,2)=exp(u03);
compressed_deriv(:,3)=-u04.*exp(u03).*(DELTA - 1);
compressed_deriv(:,4)=-exp(u03).*(DELTA - 1);
compressed_deriv(:,5)=u04.*exp(u03).*(tilrp - tilqp.*(DELTA - 1));
compressed_deriv(:,6)=exp(u03).*(tilrp - tilqp.*(DELTA - 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{2,2});
full_vals(:,4:14)=uncompressed_deriv;
full_cols(4:14,:)=index.loc{2,2};
full_rows(4:14)=2;
compressed_deriv=zeros(n_s,8);
compressed_deriv(:,1)=exp(u17 - logtilqe).*exp(logtilqep);
compressed_deriv(:,2)=-exp(u17 - logtilqe).*exp(logtilqep);
compressed_deriv(:,3)=exp(u17 - logtilqe).*exp(logtilqep);
compressed_deriv(:,4)=exp(u17 - logtilqe).*(tildivp + exp(logtilqep));
compressed_deriv(:,5)=-exp(u17 - logtilqe).*(tildivp + exp(logtilqep));
compressed_deriv(:,6)=-exp(u17 - logtilqe);
compressed_deriv(:,7)=exp(u17 - logtilqe).*(tildivp + exp(logtilqep));
compressed_deriv(:,8)=exp(u17 - logtilqe);
uncompressed_deriv=compressed_deriv(:,index.nnz{3,2});
full_vals(:,15:27)=uncompressed_deriv;
full_cols(15:27,:)=index.loc{3,2};
full_rows(15:27)=3;
compressed_deriv=zeros(n_s,3);
compressed_deriv(:,1)=exp(logMp - logqf);
compressed_deriv(:,2)=-exp(logMp - logqf);
compressed_deriv(:,3)=exp(logMp - logqf);
uncompressed_deriv=compressed_deriv(:,index.nnz{4,2});
full_vals(:,28:31)=uncompressed_deriv;
full_cols(28:31,:)=index.loc{4,2};
full_rows(28:31)=4;
derivs=sptensor(full_rows,full_cols,full_vals,4,repmat(53,1,2));
