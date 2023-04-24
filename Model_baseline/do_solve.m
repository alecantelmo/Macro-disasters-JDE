% clear,

clc
load('models')
load('do_methods')
% specify the required orders (you must have a 3rd order perturbation,
% because it is the initial guess for the other methods)

% all solutions
Orders_for_Perturbation=3;%[1,2,3,4,5];
Orders_for_TaylorProjection=[2,3];
Orders_for_Smolyak=1;%[1,2,3];

% only cheap solutions

if exist('options.mat','file')~=0
    load('options');
else
    do_all_solutions=1;
end

if do_all_solutions==0
    Orders_for_Perturbation=[1,2,3];
    Orders_for_TaylorProjection=[1,2,3];
    Orders_for_Smolyak=[1,2];
end

do_pert=do_methods(1); % do perturbation        
do_tp=do_methods(2); % do Taylor projection     
do_smol=do_methods(3); % do Smolyak collocation

main_folder=pwd;
addpath(main_folder)

% choose parameterization (1 for disaster parameters, 0 for no-disasters)
ver_options=[1]; % do disaster parameters
% ver_options=[0]; % do no-disaster parameters
% ver_options=[1,0]; % do both

for veri=ver_options
    ver=num2str(veri);

    param_folder=['params_SS' num2str(ver)];

    % Perturbation
    pert_folder=['solutions' num2str(ver) '_pert']; % folder to save results
    pert_order=Orders_for_Perturbation; % do these solution orders
    pert_algo='gensylv';

    % Taylor projection
    tp_folder=['solutions' num2str(ver) '_tp_newton']; % folder to save results
    % tols for Newton solver
    tp_tolX=1e-5;
    tp_tolF=1e-5;
    tp_maxiter=100;
    tp_jacobian='exact'; % use an exact analytic Jacobian.
    tp_order=Orders_for_TaylorProjection; % do these solution orders
    tp_maxload=50; % control the maximum data in vectorized operations. If memory is insufficient try to reduce it.

    % Smolyak
    smol_folder=['solutions' num2str(ver) '_smol_newton']; % folder to save results
    % tols for Newton solver
    smol_tolX=1e-5;
    smol_tolF=1e-5;
    smol_maxiter=100;
    smol_order=Orders_for_Smolyak; % do these solution orders
    
    grid_maxload=10; % control the maximum data in vectorized operations. If memory is insufficient try to reduce it.
    smol_maxload=100;% also controls the maximum data in vectorized operations. If memory is insufficient try to reduce it.

    save(['ver' ver],'*folder')

    solve_all_models;

end

