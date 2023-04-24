function model=prepare_collocation(f,Phi,yp,y,xp,x,eta,symparams,M,varargin)
% This code prepares codes and data that are used later to implement the
% Smolyak collocation method.
% M=Smolyak approximation level. Other arguments are like prepare_tp.
%
% © Copyright, Jesus Fernandez-Villaverde and Oren Levintal, June 13, 2016.

N=0;

currentFolder=pwd;
mkdir('files')
cd 'files'

disp('Preparing Smolyak collocation...')
% SIZE VARIABLES
n_f=length(f); n_x=length(x); n_y=length(y); n_x2=length(Phi); n_x1=n_x-n_x2;
v=[yp(:); y(:); xp(:); x(:)];
n_v=length(v);

% DEFINE pi_ and u
if isempty(varargin)
    pi_=sym([]);
    u=sym([]);
elseif isempty(varargin{1})
    warning('no substitutions assumed')
    pi_=sym([]);
    u=sym([]);
else
    pi_=varargin{1}; % pi_ is the function pi(v,u)
    pi_=pi_(:);
    u=varargin{2}; % auxiliary variables that are eventually substituted out
    u=u(:);
    if ~isequal(sort(u),unique(u))
        error('substituted variables are not uniquely determined')
    end
end
    
n_u=length(u);
model.n_u=n_u;
tilf=f; % tilf(v,u) is a function of v and u
n=find_n(pi_,u); % number of substitutions needed to eliminate u

% Identify variables that affect u
u_v=u; %u_v which is u as a function of v only u(v)
uu=eye(n_u); % matrix to store which of u affect u
uv=zeros(n_u,n_v);
for k=1:n 
    u_v=subs(u_v,u,pi_); % substitute pi_ into itself n times
    uu=uu+logical(jacobian(u_v,u)~=0); % uu is a n_u-by-n_u matrix. the ij element is zero only if ui is independent of uj through all substitutions
    uv=uv+logical(jacobian(u_v,v)~=0); % uv is a n_u-by-n_v matrix. the ij element is zero only if ui is independent of vj through all substitutions
end

% Identify stochasic and predetermined functions and variables
if n_x2>0
    stochexog=find(sum(1-logical(eta(n_x1+1:end,:)==0),2)); % stochastic exogenous variables
else
    stochexog=[];
end
preexog=1:n_x2;
preexog(stochexog)=[]; % predetermined exogenous variables

stochvars=[1:n_y,2*n_y+n_x1+stochexog']; % all stochastic variables

fv=logical(logical(jacobian(tilf,v)~=0)+logical(jacobian(tilf,u)~=0)*uv~=0); % logical Jacobian of f w.r.t v 

stochfrows=find(sum( 1-logical(fv(:,stochvars)==0),2)); % stochastic rows of f
prefrows=1:n_f;
prefrows(stochfrows)=[]; % predetermined rows of f

% variables that affect the predetermined and stochastic rows of f

prefvars=find(sum(1-logical(fv(prefrows,:)==0),1));
stochfvars=find(sum(1-logical(fv(stochfrows,:)==0),1));


% Create an m file for Phi and its derivatives w.r.t x
gen_fun_vec(Phi,symparams,x,'Phi');

% Create an m file for u
newuv=jacobian(u_v,v); % calculate uv again
stochurows=find(sum( 1-logical(newuv(:,stochvars)==0),2)); % stochastic rows of u
model.stochurows=stochurows;
preurows=1:n_u;
preurows(stochurows)=[]; % predetermined rows of u
model.preurows=preurows;
preuvars=find(sum( logical(newuv(preurows,:)~=0),1)); % variables of predetermined rows of u
stochuvars=find(sum( logical(newuv(stochurows,:)~=0),1)); % variables of stochastic rows of u
gen_fun_vec(u_v(preurows),symparams,v(preuvars),'preu');
gen_fun_vec(u_v(stochurows),symparams,v(stochuvars),'stochu');
model.preuvars=preuvars;
model.stochuvars=stochuvars;

% build z and find the Jacobian of z
z=[v;u];
n_z=n_v+n_u;
zz=[eye(n_v),zeros(n_v,n_u);uv,uu];
for i=1:n_z
    temp=zz(i,:);
    temp(temp~=0)=1:nnz(temp);
    zz(i,:)=temp;
end
zz=intarray(zz);
model.maxzz=intarray(max(zz(:)));

% build fv similar to zz
fv=double(fv);
for i=1:n_f
    temp=fv(i,:);
    temp(temp~=0)=1:nnz(temp);
    fv(i,:)=temp;
end
model.fv=fv;

tilfz=jacobian(tilf,z);

% build tilfz
tilfz=double(logical(tilfz~=0)); % only direct effects
for i=1:n_f
    temp=tilfz(i,:);
    temp(temp~=0)=1:nnz(temp);
    tilfz(i,:)=temp;
end
model.tilfz=tilfz;
model.maxtilfz=intarray(max(tilfz(:)));

% tilfz=symtilfz;
% Create an m file for tilf
pretilfzvars=find(sum(1-logical(tilfz(prefrows,:)==0),1));
model.pretilfzvars=pretilfzvars; % variables in z that affect pretilf
gen_fun_vec(tilf(prefrows),symparams,z(pretilfzvars),'pretilf');

stochtilfzvars=find(sum(1-logical(tilfz(stochfrows,:)==0),1));
model.stochtilfzvars=stochtilfzvars; % variables in z that affect stochtilf
gen_fun_vec(tilf(stochfrows),symparams,z(stochtilfzvars),'stochtilf');

% Differentiate PI and create m files

PI=[v;pi_];
prePI=PI;
prePI(n_v+stochurows)=0; % stochastic rows of PI are set to zero. Note: the first n_v rows include stochastic variables but their derivatives w.r.t are ones or zeros.
stochPI=PI;
stochPI([1:n_v,n_v+preurows])=0; % nonstochastic rows are set to zero. See the note above.

[model.prePI_ind_u]=getderivs_tensor(prePI,z,N+1,symparams,'prePI'); 
[model.stochPI_ind_u]=getderivs_tensor(stochPI,z,N+1,symparams,'stochPI');

prefuvars=find(sum(uu(logical(sum(logical(jacobian(tilf(prefrows),u)~=0),1)~=0),:),1));
stochfuvars=find(sum(uu(logical(sum(logical(jacobian(tilf(stochfrows),u)~=0),1)~=0),:),1));
model.prefuvars=prefuvars;
model.stochfuvars=stochfuvars;

prefzvars=[prefvars,n_v+prefuvars];
stochfzvars=[stochfvars,n_v+stochfuvars];
model.n_prefzvars=length(prefzvars);
model.n_stochfzvars=length(stochfzvars);

[model.stochtilf_ind_u,model.stoch_n]=gen_chainderivs_tensor(tilf(stochfrows),v(stochfvars),u(stochfuvars),pi_(stochfuvars),N+1,symparams,'stochtilf');
[model.pretilf_ind_u,model.pre_n]=gen_chainderivs_tensor(tilf(prefrows),v(prefvars),u(prefuvars),pi_(prefuvars),N+1,symparams,'pretilf');


% Store variables in struct model
model.prefrows=prefrows;
model.stochfrows=stochfrows;
model.prefvars=prefvars;
model.stochfvars=stochfvars;
model.stochexog=stochexog;
model.preexog=preexog;
model.n_stochexog=length(stochexog);
model.n_preexog=length(preexog);
model.stochfzvars=stochfzvars;
model.prefzvars=prefzvars;

model.n_f=n_f; 
model.n_x=n_x; 
model.n_y=n_y; 
model.n_x1=n_x1;
model.n_x2=n_x2; 

model.n_v=n_v;
model.n_z=n_z;
model.n_u=n_u;

model.n_ind=2;

model.stochzz=zz(model.stochfzvars,model.stochfzvars);
model.prezz=zz(model.prefzvars,model.prefzvars);

model.stochtilfz=tilfz(model.stochfrows,model.stochfzvars);
model.pretilfz=tilfz(model.prefrows,model.prefzvars);


% do smolyak using codes written by Lilia Maliar, Serguei Maliar and Rafael Valero

d=n_x;
mu=M;
model.Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu);
symx=[];
for i=1:n_x
    eval(['symx=[symx, sym(' '''' 'symx_' num2str(i) '''' ')];']); 
end

Smol_bases = sym_Smolyak_Polynomial(symx,d,mu,model.Smolyak_elem_iso);
Tcheb=Smol_bases(:);
symx=symx(:);

Tcheb=simplify(Tcheb);

nparams=length(Tcheb);
model.nparams=nparams;
model.n_theta=nparams*n_f;


Tcheb=Tcheb(:);
gen_fun_vec(Tcheb,[],symx,'Tcheb');

[model.Tcheb_ind]=getderivs_tensor(Tcheb,symx,1,[],'Tcheb'); 

disp('Smolyak collocation prepared successfully.')
rehash;
cd(currentFolder)
