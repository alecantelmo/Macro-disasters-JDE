function function_f=preu_fun(variables_v,parameters)
BETA=parameters(1);
GAMMA=parameters(2);
ALPHA=parameters(3);
RHO=parameters(4);
DELTA=parameters(5);
SIGMA=parameters(6);
logc=variables_v(1,:);
logkp=variables_v(2,:);
logk=variables_v(3,:);
loga=variables_v(4,:);
function_f=zeros(6,size(variables_v,2));
function_f(1,:)=loga + ALPHA.*logk;
function_f(2,:)=exp(loga + ALPHA.*logk);
function_f(3,:)=exp(logk);
function_f(4,:)=exp(logc);
function_f(5,:)=exp(logkp);
function_f(6,:)=exp(-logk);