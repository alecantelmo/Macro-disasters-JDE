function derivs=prePI_d4(vars,params,index)
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
full_rows=zeros(24,1);
full_cols=zeros(24,4);
full_vals=zeros(n_s,24);
compressed_deriv=zeros(n_s,2);
compressed_deriv(:,1)=exp(logtheta).*(ALPHA - 1);
compressed_deriv(:,2)=d.*exp(logtheta).*(ALPHA - 1);
uncompressed_deriv=compressed_deriv(:,index.nnz{17,4});
full_vals(:,1:5)=uncompressed_deriv;
full_cols(1:5,:)=index.loc{17,4};
full_rows(1:5)=17;
compressed_deriv=zeros(n_s,2);
compressed_deriv(:,1)=(BETA.*exp((auxvar1.*(PSI - 1))./(GAMMA - 1)).*(PSI - 1).^4)./(GAMMA - 1).^4;
compressed_deriv(:,2)=SCALEPARAM.*exp(u2);
uncompressed_deriv=compressed_deriv(:,index.nnz{21,4});
full_vals(:,6:7)=uncompressed_deriv;
full_cols(6:7,:)=index.loc{21,4};
full_rows(6:7)=21;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=6./(u1.^4.*(PSI - 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{25,4});
full_vals(:,8:8)=uncompressed_deriv;
full_cols(8:8,:)=index.loc{25,4};
full_rows(8:8)=25;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(nl - log(exp(nl) + 1)).*(exp(nl)./(exp(nl) + 1) - 1).^4 + exp(nl - log(exp(nl) + 1)).*((7.*exp(2.*nl))./(exp(nl) + 1).^2 - (12.*exp(3.*nl))./(exp(nl) + 1).^3 + (6.*exp(4.*nl))./(exp(nl) + 1).^4 - exp(nl)./(exp(nl) + 1)) + 3.*exp(nl - log(exp(nl) + 1)).*(exp(2.*nl)./(exp(nl) + 1).^2 - exp(nl)./(exp(nl) + 1)).^2 + 4.*exp(nl - log(exp(nl) + 1)).*(exp(nl)./(exp(nl) + 1) - 1).*((2.*exp(3.*nl))./(exp(nl) + 1).^3 - (3.*exp(2.*nl))./(exp(nl) + 1).^2 + exp(nl)./(exp(nl) + 1)) + 6.*exp(nl - log(exp(nl) + 1)).*(exp(nl)./(exp(nl) + 1) - 1).^2.*(exp(2.*nl)./(exp(nl) + 1).^2 - exp(nl)./(exp(nl) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{27,4});
full_vals(:,9:9)=uncompressed_deriv;
full_cols(9:9,:)=index.loc{27,4};
full_rows(9:9)=27;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-(6.*NU)./(l - 1).^4;
uncompressed_deriv=compressed_deriv(:,index.nnz{29,4});
full_vals(:,10:10)=uncompressed_deriv;
full_cols(10:10,:)=index.loc{29,4};
full_rows(10:10)=29;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-(6.*NU)./(l - 1).^4;
uncompressed_deriv=compressed_deriv(:,index.nnz{31,4});
full_vals(:,11:11)=uncompressed_deriv;
full_cols(11:11,:)=index.loc{31,4};
full_rows(11:11)=31;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-6./(l - 1).^4;
uncompressed_deriv=compressed_deriv(:,index.nnz{35,4});
full_vals(:,12:12)=uncompressed_deriv;
full_cols(12:12,:)=index.loc{35,4};
full_rows(12:12)=35;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logtheta);
uncompressed_deriv=compressed_deriv(:,index.nnz{40,4});
full_vals(:,13:13)=uncompressed_deriv;
full_cols(13:13,:)=index.loc{40,4};
full_rows(13:13)=40;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-(ALPHA - 1).*((7.*exp(2.*nl))./(exp(nl) + 1).^2 - (12.*exp(3.*nl))./(exp(nl) + 1).^3 + (6.*exp(4.*nl))./(exp(nl) + 1).^4 - exp(nl)./(exp(nl) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{45,4});
full_vals(:,14:14)=uncompressed_deriv;
full_cols(14:14,:)=index.loc{45,4};
full_rows(14:14)=45;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-ALPHA.*((7.*exp(2.*nl))./(exp(nl) + 1).^2 - (12.*exp(3.*nl))./(exp(nl) + 1).^3 + (6.*exp(4.*nl))./(exp(nl) + 1).^4 - exp(nl)./(exp(nl) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{47,4});
full_vals(:,15:15)=uncompressed_deriv;
full_cols(15:15,:)=index.loc{47,4};
full_rows(15:15)=47;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-(ALPHA - 1).*((7.*exp(2.*nl))./(exp(nl) + 1).^2 - (12.*exp(3.*nl))./(exp(nl) + 1).^3 + (6.*exp(4.*nl))./(exp(nl) + 1).^4 - exp(nl)./(exp(nl) + 1));
uncompressed_deriv=compressed_deriv(:,index.nnz{49,4});
full_vals(:,16:16)=uncompressed_deriv;
full_cols(16:16,:)=index.loc{49,4};
full_rows(16:16)=49;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logtilc);
uncompressed_deriv=compressed_deriv(:,index.nnz{50,4});
full_vals(:,17:17)=uncompressed_deriv;
full_cols(17:17,:)=index.loc{50,4};
full_rows(17:17)=50;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(u4);
uncompressed_deriv=compressed_deriv(:,index.nnz{52,4});
full_vals(:,18:18)=uncompressed_deriv;
full_cols(18:18,:)=index.loc{52,4};
full_rows(18:18)=52;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-6./tilx.^4;
uncompressed_deriv=compressed_deriv(:,index.nnz{54,4});
full_vals(:,19:19)=uncompressed_deriv;
full_cols(19:19,:)=index.loc{54,4};
full_rows(19:19)=54;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(-auxvar1);
uncompressed_deriv=compressed_deriv(:,index.nnz{57,4});
full_vals(:,20:20)=uncompressed_deriv;
full_cols(20:20,:)=index.loc{57,4};
full_rows(20:20)=57;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=24./tilq.^5;
uncompressed_deriv=compressed_deriv(:,index.nnz{59,4});
full_vals(:,21:21)=uncompressed_deriv;
full_cols(21:21,:)=index.loc{59,4};
full_rows(21:21)=59;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(-logtilx);
uncompressed_deriv=compressed_deriv(:,index.nnz{60,4});
full_vals(:,22:22)=uncompressed_deriv;
full_cols(22:22,:)=index.loc{60,4};
full_rows(22:22)=60;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=-6./tily.^4;
uncompressed_deriv=compressed_deriv(:,index.nnz{65,4});
full_vals(:,23:23)=uncompressed_deriv;
full_cols(23:23,:)=index.loc{65,4};
full_rows(23:23)=65;
compressed_deriv=zeros(n_s,1);
compressed_deriv(:,1)=exp(logtilw);
uncompressed_deriv=compressed_deriv(:,index.nnz{68,4});
full_vals(:,24:24)=uncompressed_deriv;
full_cols(24:24,:)=index.loc{68,4};
full_rows(24:24)=68;
derivs=sptensor(full_rows,full_cols,full_vals,71,repmat(71,1,4));