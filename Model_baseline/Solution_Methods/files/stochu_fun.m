function function_f=stochu_fun(variables_v,parameters)
BETA=parameters(1);
GAMMA=parameters(2);
ALPHA=parameters(3);
RHO=parameters(4);
DELTA=parameters(5);
SIGMA=parameters(6);
logcp=variables_v(1,:);
logc=variables_v(2,:);
logkp=variables_v(3,:);
logap=variables_v(4,:);
function_f=zeros(4,size(variables_v,2));
function_f(1,:)=exp(logap + log(ALPHA) + logkp.*(ALPHA - 1));
function_f(2,:)=logap + log(ALPHA) + logkp.*(ALPHA - 1);
function_f(3,:)=exp(log(BETA) + GAMMA.*(logc - logcp));
function_f(4,:)=log(BETA) + GAMMA.*(logc - logcp);
