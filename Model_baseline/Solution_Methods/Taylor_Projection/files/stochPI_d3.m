function derivs=stochPI_d3(vars,params,index)
n_s=size(vars,1);
auxvar1p=vars(:,1);
nlp=vars(:,2);
logtilqep=vars(:,3);
logqfp=vars(:,4);
auxvar1=vars(:,5);
nl=vars(:,6);
logtilqe=vars(:,7);
logqf=vars(:,8);
logtilkstarbackp=vars(:,9);
dp=vars(:,10);
logthetap=vars(:,11);
zap=vars(:,12);
logtilkstarback=vars(:,13);
d=vars(:,14);
logtheta=vars(:,15);
za=vars(:,16);
loghata=vars(:,17);
loghatap=vars(:,18);
loghatz=vars(:,19);
loghatzp=vars(:,20);
u1=vars(:,21);
u1p=vars(:,22);
u2=vars(:,23);
u2p=vars(:,24);
logtilv_tilvss=vars(:,25);
logtilvp_tilvss=vars(:,26);
l=vars(:,27);
lp=vars(:,28);
logtilu=vars(:,29);
logtilup=vars(:,30);
loguc=vars(:,31);
logucp=vars(:,32);
logtillambda=vars(:,33);
logtillambdap=vars(:,34);
logtilc=vars(:,35);
logtilcp=vars(:,36);
logmnom=vars(:,37);
logmnomp=vars(:,38);
logmden=vars(:,39);
theta=vars(:,40);
thetap=vars(:,41);
tilrp=vars(:,42);
tilq=vars(:,43);
tilqp=vars(:,44);
logtilr=vars(:,45);
logtilrp=vars(:,46);
logtilw=vars(:,47);
logtilwp=vars(:,48);
u4=vars(:,49);
tilc=vars(:,50);
tilcp=vars(:,51);
tilx=vars(:,52);
tilxp=vars(:,53);
logtilx=vars(:,54);
logtilk=vars(:,55);
u01p=vars(:,56);
u02=vars(:,57);
u03=vars(:,58);
u04=vars(:,59);
u05=vars(:,60);
logMp=vars(:,61);
u17=vars(:,62);
u4p=vars(:,63);
tily=vars(:,64);
logtily=vars(:,65);
tilyp=vars(:,66);
logtilyp=vars(:,67);
w=vars(:,68);
wp=vars(:,69);
tildiv=vars(:,70);
tildivp=vars(:,71);
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
full_rows=zeros(21,1);
full_cols=zeros(21,3);
full_vals=zeros(n_s,21);
compressed_deriv=zeros(n_s,2);
compressed_deriv(:,1)=exp(logthetap).*(ALPHA - 1);
compressed_deriv(:,2)=dp.*exp(logthetap).*(ALPHA - 1);
uncompressed_deriv=compressed_deriv(:,index.nnz{18,3});
full_vals(:,1:4)=uncompressed_deriv;
full_cols(1:4,:)=index.loc{18,3};
full_rows(1:4)=18;
compressed_deriv=zeros(n_s,2);
compressed_deriv(:,1)=(BETA.*exp((auxvar1p.*(PSI - 1))./(GAMMA - 1)).*(PSI - 1).^3)./(GAMMA - 1).^3;
compressed_deriv(:,2)=SCALEPARAM.*exp(u2p);
uncompressed_deriv=compressed_deriv(:,index.nnz{22,3});
full_vals(:,5:6)=uncompressed_deriv;
full_cols(5:6,:)=index.loc{22,3};
full_rows(5:6)=22;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-2./(u1p.^3.*(PSI - 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{26,3});
full_vals(:,7:7)=uncompressed_deriv;
full_cols(7:7,:)=index.loc{26,3};
full_rows(7:7)=26;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=- exp(nlp - log(exp(nlp) + 1)).*((2.*exp(3.*nlp))./(exp(nlp) + 1).^3 - (3.*exp(2.*nlp))./(exp(nlp) + 1).^2 + exp(nlp)./(exp(nlp) + 1)) - exp(nlp - log(exp(nlp) + 1)).*(exp(nlp)./(exp(nlp) + 1) - 1).^3 - 3.*exp(nlp - log(exp(nlp) + 1)).*(exp(nlp)./(exp(nlp) + 1) - 1).*(exp(2.*nlp)./(exp(nlp) + 1).^2 - exp(nlp)./(exp(nlp) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{28,3});
full_vals(:,8:8)=uncompressed_deriv;
full_cols(8:8,:)=index.loc{28,3};
full_rows(8:8)=28;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=(2.*NU)./(lp - 1).^3;
uncompressed_deriv=compressed_deriv(:,index.nnz{30,3});
full_vals(:,9:9)=uncompressed_deriv;
full_cols(9:9,:)=index.loc{30,3};
full_rows(9:9)=30;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=(2.*NU)./(lp - 1).^3;
uncompressed_deriv=compressed_deriv(:,index.nnz{32,3});
full_vals(:,10:10)=uncompressed_deriv;
full_cols(10:10,:)=index.loc{32,3};
full_rows(10:10)=32;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=2./(lp - 1).^3;
uncompressed_deriv=compressed_deriv(:,index.nnz{36,3});
full_vals(:,11:11)=uncompressed_deriv;
full_cols(11:11,:)=index.loc{36,3};
full_rows(11:11)=36;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logthetap);
uncompressed_deriv=compressed_deriv(:,index.nnz{41,3});
full_vals(:,12:12)=uncompressed_deriv;
full_cols(12:12,:)=index.loc{41,3};
full_rows(12:12)=41;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logtilrp);
uncompressed_deriv=compressed_deriv(:,index.nnz{42,3});
full_vals(:,13:13)=uncompressed_deriv;
full_cols(13:13,:)=index.loc{42,3};
full_rows(13:13)=42;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=(ALPHA - 1).*((2.*exp(3.*nlp))./(exp(nlp) + 1).^3 - (3.*exp(2.*nlp))./(exp(nlp) + 1).^2 + exp(nlp)./(exp(nlp) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{46,3});
full_vals(:,14:14)=uncompressed_deriv;
full_cols(14:14,:)=index.loc{46,3};
full_rows(14:14)=46;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=ALPHA.*((2.*exp(3.*nlp))./(exp(nlp) + 1).^3 - (3.*exp(2.*nlp))./(exp(nlp) + 1).^2 + exp(nlp)./(exp(nlp) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{48,3});
full_vals(:,15:15)=uncompressed_deriv;
full_cols(15:15,:)=index.loc{48,3};
full_rows(15:15)=48;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logtilcp);
uncompressed_deriv=compressed_deriv(:,index.nnz{51,3});
full_vals(:,16:16)=uncompressed_deriv;
full_cols(16:16,:)=index.loc{51,3};
full_rows(16:16)=51;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(u4p);
uncompressed_deriv=compressed_deriv(:,index.nnz{53,3});
full_vals(:,17:17)=uncompressed_deriv;
full_cols(17:17,:)=index.loc{53,3};
full_rows(17:17)=53;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=(ALPHA - 1).*((2.*exp(3.*nlp))./(exp(nlp) + 1).^3 - (3.*exp(2.*nlp))./(exp(nlp) + 1).^2 + exp(nlp)./(exp(nlp) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{63,3});
full_vals(:,18:18)=uncompressed_deriv;
full_cols(18:18,:)=index.loc{63,3};
full_rows(18:18)=63;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(u4p);
uncompressed_deriv=compressed_deriv(:,index.nnz{66,3});
full_vals(:,19:19)=uncompressed_deriv;
full_cols(19:19,:)=index.loc{66,3};
full_rows(19:19)=66;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=2./tilyp.^3;
uncompressed_deriv=compressed_deriv(:,index.nnz{67,3});
full_vals(:,20:20)=uncompressed_deriv;
full_cols(20:20,:)=index.loc{67,3};
full_rows(20:20)=67;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logtilwp);
uncompressed_deriv=compressed_deriv(:,index.nnz{69,3});
full_vals(:,21:21)=uncompressed_deriv;
full_cols(21:21,:)=index.loc{69,3};
full_rows(21:21)=69;
derivs=sptensor(full_rows,full_cols,full_vals,71,repmat(71,1,3));
