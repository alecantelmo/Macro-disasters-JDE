function function_f=pretilf_fun(variables_v,parameters)
MUD=parameters(1);
RHOTHETA=parameters(2);
THETABAR=parameters(3);
SDV_THETA=parameters(4);
SDV_ZA=parameters(5);
ALPHA=parameters(6);
LAMBDA_A=parameters(7);
GAMMA=parameters(8);
NU=parameters(9);
PSI=parameters(10);
BETA=parameters(11);
DELTA=parameters(12);
PHI=parameters(13);
SCALEPARAM=parameters(14);
LOGTILUSS=parameters(15);
logtilkstarbackp=variables_v(1,:);
logtilk=variables_v(2,:);
u05=variables_v(3,:);
function_f=zeros(1,size(variables_v,2));
function_f(1,:)=u05.*exp(logtilkstarbackp) + u05.*exp(logtilk).*(DELTA - 1) - 1;
