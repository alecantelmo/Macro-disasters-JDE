function derivs=Phi_d3(vars,params,index)
logtilkstarback=vars(1);
d=vars(2);
logtheta=vars(3);
za=vars(4);
sigma_perturbation=vars(5);
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
full_rows=zeros(0,1);
full_cols=zeros(0,1);
full_vals=zeros(0,1);
derivs=sparse(full_rows,full_cols,full_vals,3,125);