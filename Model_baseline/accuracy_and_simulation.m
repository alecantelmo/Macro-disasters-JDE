% This file performs the following:
% 1. Computes model residuals over the ergodic set. The ergodic set is
%    approximated by simulating the model for 1000 periods with the 3rd order TP
%    solution.
% 2. Simulates the model from the steady state using all solutions.

% clear,
clc, close all
load('models') 
load('do_methods')
% models=models(7); %take a specific model

%%% all solutions
Orders_for_Perturbation=3;%[1,2,3,4,5];
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

do_accuracy=1; % compute model residuals over the ergodic set
do_simulation=1; % simulate the model by the different solutions

main_folder=pwd;
addpath(main_folder)

all_versions={'ver1'}; % do disaster parameterizaion
% all_versions={'ver0'}; % do no-disaster parameterization
% all_versions={'ver1','ver0'}; % do both

for veri=1:length(all_versions)
    version=all_versions{veri};
    load(version)

    for modeli=1:length(models)
        mkdir(['performance\' version '\' models{modeli}])

        % Take the model equations for accuracy tests from this path:
        addpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\files']);
        load([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\model'],'model');
        load([main_folder '\all_models\' models{modeli} '\TaylorProjection\choose']);

        n_y=model.n_y;
        n_x=model.n_x;
        n_f=model.n_f;

        % Approximate the ergodic set by tp3
        load([main_folder '\' tp_folder '\' models{modeli} '\order3\solution_tp']);
        rng('default');
        T=1000;
        load([main_folder '\' pert_folder '\' models{modeli} '\order3\perturbation'],'eta_mat','nyss','nxss','params','P1','nep1','n_e','prob_disaster','MUD')

        shocks=[(rand(1,T)<prob_disaster)-MUD;randn(n_e-1,T)];

        N=3; eta=eta_mat; [U3,W3]=create_UW(n_x,3); model.U{3}=U3;model.W{3}=W3;
        [ results,nu ] = fun_simulation( nxss,shocks,coeffs,c0,N,n_y,n_x,n_f,model,params,eta );
        
        % this is the ergodic set
        xt_simul=results(n_y+1:end,:);
        
        % load the basic model specifications
        load([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\model'],'model');
        model0=model;

        n_y=model0.n_y;
        n_x1=model0.n_x1;
        n_x=model0.n_x;
        n_f=model0.n_f;

        % load Taylor Projection
        for order=Orders_for_TaylorProjection
            load([main_folder '\' tp_folder '\' models{modeli} '\order' num2str(order) '\solution_tp'],'coeffs','time');
            coeffs=reshape(coeffs,model0.n_f,numel(coeffs)/model0.n_f); coeffs=coeffs(:);
            eval(['tp' num2str(order) '=coeffs;']) 
            eval(['tp' num2str(order) '_time=time;'])
        end

        % load Smolyak
     if do_methods(3)==1
        for order=Orders_for_Smolyak
            load([main_folder '\' smol_folder '\' models{modeli} '\order' num2str(order) '\solution_smol'],'coeffs','cheb_c0','DELTA','time');
            coeffs=reshape(coeffs,model0.n_f,numel(coeffs)/(model0.n_f)); coeffs=coeffs(:);
            eval(['smol' num2str(order) '=coeffs;'])
            eval(['smol' num2str(order) '_c0=cheb_c0;'])
            eval(['smol' num2str(order) '_DELTA=DELTA;'])
            eval(['smol' num2str(order) '_n_f=model.n_f;'])
            eval(['smol' num2str(order) '_nparams=model.n_b;'])
            eval(['smol' num2str(order) '_time=time;'])
        end
     end
        % load perturbation runtime
        for order=Orders_for_Perturbation
            load([main_folder '\' pert_folder '\' models{modeli} '\order' num2str(order) '\perturbation'],'time');
            eval(['pert' num2str(order) '_time=time;'])
        end
    
        % create perturbation solutions of orders 1,2,3,4,5 by changing the
        % perturbation parameters from 0 to 1.
        
        for order=Orders_for_Perturbation
            % load perturbation solution
            load([main_folder '\' pert_folder '\' models{modeli} '\order' num2str(order) '\perturbation'],'derivs','nxss','nyss','params','eta_mat','nep_full','P_full','P1','nep1');

            nyss=nyss(choose);
            eta=eta_mat;

            pert_coeffs = GetInitial( order,derivs,choose,nyss,nxss,n_x,n_x1 );
            eval(['pert' num2str(order) '=pert_coeffs;']);
        end
        
        T=size(shocks,2);
        nparams=model.n_b;
        tp1_r=zeros(n_f,T);
        tp2_r=zeros(n_f,T);
        tp3_r=zeros(n_f,T);

        pert1_r=zeros(n_f,T);
        pert2_r=zeros(n_f,T);
        pert3_r=zeros(n_f,T);
        pert4_r=zeros(n_f,T);
        pert5_r=zeros(n_f,T);

        smol1_r=zeros(n_f,T);
        smol2_r=zeros(n_f,T);
        smol3_r=zeros(n_f,T);

        et=shocks;
        c0=nxss;
        model=model0;
        
        [~,W2]=create_UW(n_x,2);
        model.W{2}=W2;
        [~,W3]=create_UW(n_x,3);
        model.W{3}=W3;
        [~,W4]=create_UW(n_x,4);
        model.W{4}=W4;
        [~,W5]=create_UW(n_x,5);
        model.W{5}=W5;

        
        if do_accuracy==1
            P=P1;
            nep=nep1;
            for t=1:T
                t
                x0=xt_simul(:,t);

                
                % Taylor Projection
                poly.type='power';
                for order=Orders_for_TaylorProjection
                    eval(['coeffs=tp' num2str(order) '(:);'])
                    resids=eval_eqs(coeffs,x0,model,params,eta,c0,nep,P,order,poly);
                    eval(['tp' num2str(order) '_r(:,t)=resids;'])
                end
                
                % Perturbation
                poly.type='power';
                for order=Orders_for_Perturbation
                    eval(['coeffs=pert' num2str(order) '(:);'])
                    resids=eval_eqs(coeffs,x0,model,params,eta,c0,nep,P,order,poly);
                    eval(['pert' num2str(order) '_r(:,t)=resids;']);
                end
                
                % smolyak
                if do_methods(3)==1 %
                for order=Orders_for_Smolyak
                    poly.type='smol';
                    eval(['poly.DELTA=smol' num2str(order) '_DELTA;'])
                    poly.Tcheb_fun_folder=[main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order) '\files'];
                    eval(['poly.n_f=smol' num2str(order) '_n_f;'])
                    eval(['poly.nparams=smol' num2str(order) '_nparams;'])
                    eval(['coeffs=smol' num2str(order) '(:);'])
                    eval(['smol_c0=smol' num2str(order) '_c0;'])
                    resids=eval_eqs(coeffs,x0,model,params,eta,smol_c0,nep,P,order,poly);
                    eval(['smol' num2str(order) '_r(:,t)=resids;']);
                end
                end
            end
                save([main_folder '\performance\' version '\' models{modeli} '\accuracy_and_runtime'],'pert*','tp*','smol*','shocks','T')
            else
                load([main_folder '\performance\' version '\' models{modeli} '\accuracy_and_runtime'])
        end


        takevars=[1:n_y];
        

        if do_simulation==1
            % Taylor projection
            for order=Orders_for_TaylorProjection
                eval(['coeffs=tp' num2str(order) '(:);'])
                [ results,nu ] = fun_simulation( nxss,shocks,coeffs,c0,order,n_y,n_x,n_f,model,params,eta );
                eval(['results_tp' num2str(order) '=results(takevars,:);'])
                eval(['x_tp' num2str(order) '=results(n_y+1:end,:);'])
                eval(['nu_tp' num2str(order) '=nu;'])
            end

            % perturbation
            for order=Orders_for_Perturbation
                eval(['coeffs=pert' num2str(order) '(:);'])
                [ results,nu] = fun_simulation( nxss,shocks,coeffs,c0,order,n_y,n_x,n_f,model,params,eta );
                eval(['results_pert' num2str(order) '=results(takevars,:);'])
                eval(['x_pert' num2str(order) '=results(n_y+1:end,:);'])
                eval(['nu_pert' num2str(order) '=nu;'])
            end
            
            % Smolyak
            if do_methods(3)==1 %added by Alessandro
            for order=Orders_for_Smolyak
                eval(['coeffs=smol' num2str(order) '(:);'])
                eval(['smol_c0=smol' num2str(order) '_c0;'])
                poly.type='smol';
                eval(['poly.DELTA=smol' num2str(order) '_DELTA;'])

                poly.Tcheb_fun_folder=[main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order) '\files'];
                [ results,nu] = fun_simulation( nxss,shocks,coeffs,smol_c0,order,n_y,n_x,n_f,model,params,eta,poly );
                eval(['results_smol' num2str(order) '=results(takevars,:);'])
                eval(['x_smol' num2str(order) '=results(n_y+1:end,:);'])
                eval(['nu_smol' num2str(order) '=nu;'])
            end
            else
            x_smol=0;
       
            end
        

            save([main_folder '\performance\' version '\' models{modeli} '\simulation'],'results_*','x_tp*','x_pert*','nu_*')
        else
            load([main_folder '\performance\' version '\' models{modeli} '\simulation'],'results_*','x_tp*','x_pert*','x_smol*','nu_*')
        end
        rmpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\files']);

    end
end

