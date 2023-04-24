% This file prepares codes and data that are used later to solve the eight
% models by perturbation, Taylor projection and Smolyak collocation.
clear
clc

models={%'RBC_EZ';
%     'RBC_EZ_adjCost';
    'RBC_EZ_adjCost_Calvo'};
%     'RBC_EZ_adjCost_Calvo_Taylor1';
%     'RBC_EZ_adjCost_Calvo_Taylor2';
%     'RBC_EZ_adjCost_Calvo_Taylor2_shock1';
%     'RBC_EZ_adjCost_Calvo_Taylor2_shock2';
%     'RBC_EZ_adjCost_Calvo_Taylor2_shock3'};

save('models','models')

% specify the required orders (you must have a 3rd order perturbation and a
% 2nd order Taylor projection, because the other algorithms assume these
% solutions exist).

%%% all solutions
Orders_for_Perturbation=[3,5];%[1,2,3,4,5];
Orders_for_TaylorProjection=[2,3];
Orders_for_Smolyak=1;%[1,2,3];

%%% only cheap solutions

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

do_pert=1; % prepare perturbation
do_tp=1; % prepare Taylor projection
do_smol=0; % prepare Smolyak collocation

do_methods= [do_pert, do_tp, do_smol]; 
save('do_methods','do_methods');       

main_folder=pwd;
addpath(main_folder)

% prepare perturbation

if do_pert==1
    mkdir([main_folder '\files_for_perturbation']);
    for modeli=1:length(models)
        addpath([main_folder '\all_models\' models{modeli} '\Perturbation']);
        get_model_pert;
        mkdir([main_folder '\files_for_perturbation\' models{modeli}]);
        for order=Orders_for_Perturbation
            mkdir([main_folder '\files_for_perturbation\' models{modeli} '\order' num2str(order)]);
            cd([main_folder '\files_for_perturbation\' models{modeli} '\order' num2str(order)]);
            model=differentiate_dsge(f,yp,y,xp,x,symparams,order,Phi_fun);
            save('model')
        end
        rmpath([main_folder '\all_models\' models{modeli} '\Perturbation']);
    end
end

cd(main_folder)

% prepare Taylor projection

if do_tp==1
    mkdir([main_folder '\files_for_TaylorProjection']);
    for modeli=1:length(models)
        addpath([main_folder '\all_models\' models{modeli} '\TaylorProjection']);
        get_model_tp;
        mkdir([main_folder '\files_for_TaylorProjection\' models{modeli}]);
        for order=Orders_for_TaylorProjection
            % create a new folder
            mkdir([main_folder '\files_for_TaylorProjection\' models{modeli} '\order' num2str(order)]);
            cd([main_folder '\files_for_TaylorProjection\' models{modeli} '\order' num2str(order)]);
            % prepare and save model
            model=prepare_tp(f,Phi_fun,yp,y,xp,x,symparams,eta_mat,order,subsfuns,subsvars);
            save('model')
        end
        rmpath([main_folder '\all_models\' models{modeli} '\TaylorProjection']);
    end
end

cd(main_folder)


% prepare Smolyak collocation

if do_smol==1
    mkdir([main_folder '\files_for_smolyak']);
    for modeli=1:length(models)
        addpath([main_folder '\all_models\' models{modeli} '\TaylorProjection']);
        get_model_tp; % the model is the same as the model used for Taylor Projection
        mkdir([main_folder '\files_for_smolyak\' models{modeli}]);
        for order=Orders_for_Smolyak
            mkdir([main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order)]);
            cd([main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order)]);
            % prepare and save model
            model=prepare_collocation(f,Phi_fun,yp,y,xp,x,eta_mat,symparams,order,subsfuns,subsvars);
            save('model')
        end
        rmpath([main_folder '\all_models\' models{modeli} '\TaylorProjection']);
    end
end

cd(main_folder)