function model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi)
%differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi,ufun,u)
%This function differentiates the dsge model and prepares all the files 
%necessary to run solve_dsge.m.
%Input variables:
%f,yp,y,xp,x: the model conditions and variable as in Schmitt-Grohe and Uribe (2004)
%symparams: a symbolic array that lists all parameters.
%approx: order of perturbation
%Phi: expected value of exgoenous variables (leave empty if not
%used)
%
% © Copyright, Oren Levintal, June 13, 2016.


disp('Preparing perturbation...')

currentFolder = pwd;
mkdir('files');
cd('files')

syms sigma_perturbation sigmap_perturbation real
Phi=sym(Phi);
f=sym(f);

x=[x(:);sigma_perturbation]; 
xp=[xp(:);sigmap_perturbation];

v=[yp(:); y(:); xp(:); x(:)];

n_f=length(f); 
n_x=length(x); 
n_y=length(y); 
n_x2=size(Phi,1); 
n_x1=n_f-n_y;
n_v=length(v);

model.f=f;
model.yp=yp;
model.y=y;
model.xp=xp;
model.x=x;
model.v=v;
model.n_f=n_f;
model.n_x=n_x;
model.n_y=n_y;
model.n_x2=n_x2;
model.n_x1=n_x1;
model.n_v=n_v;

create_compression_matrices_nonzero;

if approx>=1
    model.f_ind=getderivs_vector(f,v,approx,symparams,'f');
end

% Create an m file for Phi and its derivatives w.r.t x

gen_fun(Phi,symparams,x,'Phi');
if approx>=1
    model.Phi_ind=getderivs(Phi,x,approx,symparams,'Phi');
end

create_OMEGA_x;

model.UW=UW;

if approx<3
    OMEGA_x=[];
end
model.OMEGA_x=OMEGA_x;

rehash;

disp('Perturbation prepared successfully.')

cd(currentFolder)

