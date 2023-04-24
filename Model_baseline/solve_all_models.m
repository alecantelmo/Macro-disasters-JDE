% Solve by Perturbation

if do_pert==1
    mkdir([main_folder '\' pert_folder]);
    for modeli=1:length(models)
        cd([main_folder '\' param_folder '\' models{modeli}]);

        % Choose parameter values 
        parameters;             %calls the file parameters in params_SS1 (after also SS0) in the corresponding model folder 
        
        % Discretize shocks
        discrete_shocks;

        P=P1; nep=nep1; % use first order Monomials to discretize Gaussian Shocks.
        
        % Define the cross moments in a structure with field names M2, M3, M4, M5
        moments_discrete; 

        % Steady state values
        fprintf(['\n' models{modeli} ': solving steady state... \n'])

        tic
        SteadyState;
        ss_time=toc;
       
        % build nyss and nxss, which are vectors of the steady state values
        % of the control and state variables, and params is a vector of the
        % parameter values.
        ss_and_params;
        
        cd(main_folder)
        for order=pert_order
            addpath([main_folder '\files_for_perturbation\' models{modeli} '\order' num2str(order) '\files']);

            load([main_folder '\files_for_perturbation\' models{modeli} '\order' num2str(order) '\model'],'model','f','y*','x*','symparams');

            % Test steady state

            if norm(double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss;nyss;nxss;nxss;params(:)])))>1e-6
                error('Steady State solution is not accurate')
            end

            % Solve
            %algo='dlyap'; 
            %algo='vectorize';
            %algo='gensylv';
            algo=pert_algo;
            

            approx=order;
            % solve first time for compilation
            fprintf(['\n' models{modeli} ': solving perturbation, order ' num2str(approx) ', first run ... \n'])

            solve_dsge(model,params,M,eta_mat,nxss,nyss,approx,algo);
            
            % solve again to measure runtime
            fprintf(['\n' models{modeli} ': solving perturbation, order ' num2str(approx) ', second run ... \n'])
            
            start=tic;
            derivs=solve_dsge(model,params,M,eta_mat,nxss,nyss,approx,algo);
            time=toc(start);

            mkdir([main_folder '\' pert_folder '\' models{modeli} '\order' num2str(order)])
            save([main_folder '\' pert_folder '\' models{modeli} '\order' num2str(order)  '\perturbation'],'derivs','model','params','eta_mat','nxss','nyss','P1','nep1','n_e','prob_disaster','MUD','time' );
            rmpath([main_folder '\files_for_perturbation\' models{modeli} '\order' num2str(order) '\files']);
        end
    end
%     rmpath(pathname);
end

cd(main_folder)

% Solve by Taylor Projection

if do_tp==1
    
    mkdir([main_folder '\' tp_folder]);
    tolX=tp_tolX;tolF=tp_tolF;maxiter=tp_maxiter; jacobian_type=tp_jacobian;
    
    for modeli=1:length(models)
%       use 3rd order perturbation for initial guess and parameters
        load([main_folder '\' pert_folder '\' models{modeli} '\order3\perturbation'],'eta_mat','params','P1','nep1','derivs','nyss','nxss')

        cd([main_folder '\all_models\' models{modeli} '\TaylorProjection']);
        do_choose;
        cd(main_folder);
        eta=eta_mat;
        nyss=nyss(choose);
        P=P1;nep=nep1;
        for order=tp_order
            addpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order' num2str(order) '\files']);
            load([main_folder '\files_for_TaylorProjection\' models{modeli} '\order' num2str(order) '\model'],'model')

            n_x=model.n_x;
            n_x1=model.n_x1;
            n_x2=model.n_x2;
            n_y=model.n_y; 
            n_f=model.n_f;

            % Make Initial Guess - third order perturbation solution
            
            [ pert_coeffs ] = GetInitial( order,derivs,choose,nyss,nxss,n_x,n_x1 );


            coeffs=pert_coeffs(:);
            c0=nxss;
            x0=nxss;
            N=order;

            % Run first time to compile files and update model
            model.maxload=tp_maxload;
            model.jacobian=jacobian_type;

            [~,~,model]=tp(coeffs,x0,model,params,eta,c0,nep,P);

            % Solve the model and save
            mkdir([main_folder '\' tp_folder '\' models{modeli} '\order' num2str(order)])
            try
                fprintf(['\n' models{modeli} ': solving Taylor projection, order ' num2str(order) '\n'])
                start=tic;
                [ncoeffs,R,J,iter]=tp_newton(coeffs,x0,model,params,eta,c0,nep,P,order,tolX,tolF,maxiter);
                time=toc(start);
                coeffs=ncoeffs;
            catch % newton failed
                ncoeffs=NaN*zeros(size(coeffs));
                R=[];
                J=[];
                iter=NaN;
                time=NaN;
            end
            save([main_folder '\' tp_folder '\' models{modeli} '\order' num2str(order)  '\solution_tp'],'coeffs','c0','R','J','iter','time');
           
            rmpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order' num2str(order) '\files']);
        end
    end
end

cd(main_folder)

% Solve by Smolyak Collocation

if do_smol==1

    mkdir([main_folder '\' smol_folder]);
    tolX=smol_tolX;tolF=smol_tolF;maxiter=smol_maxiter;
    %%%%%%%
    
    for modeli=1:length(models)
        load([main_folder '\all_models\' models{modeli} '\TaylorProjection\choose'],'choose')
        load([main_folder '\' pert_folder '\' models{modeli} '\order3\perturbation'],'eta_mat','nyss','nxss','derivs','params','P1','nep1','n_e','prob_disaster','MUD')
        eta=eta_mat;
        nyss=nyss(choose);
        
        % get the basic model
        load([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\model'],'model')
        n_x1=model.n_x1;
        n_x=model.n_x;
        % get 3rd order perturbation solution
        [ pert_coeffs ] = GetInitial( 3,derivs,choose,nyss,nxss,n_x,n_x1 );
        c0=nxss;
        simul_coeffs=pert_coeffs;
        
        % Use the perturbation solution to simulate the model in order to
        % get the approximate ergodic set
        
        addpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\files'])
        
        rng('default');
        T=1000;
        shocks=[(rand(1,T)<prob_disaster)-MUD;randn(n_e-1,T)];
        [ results ] = fun_simulation( nxss,shocks,simul_coeffs,c0,3,model.n_y,model.n_x,model.n_f,model,params,eta_mat );

        rmpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order2\files'])

        xt_simul=results(model.n_y+1:end,:);
        
        % bounds of ergodic set, mid point and DELTA for the grid
        bounds=[min(xt_simul,[],2),max(xt_simul,[],2)]; % bounds of hypercube
        cheb_c0=(bounds(:,1)+bounds(:,2))/2; % center of hypercube
        DELTA=(bounds(:,2)-bounds(:,1))/2; % radius of hypercube
        DELTA_original=DELTA;
        P=P1;nep=nep1;
        for order=smol_order
            DELTA=DELTA_original;
            % increase DELTA for 3rd order Smolyak
            if order==3
                DELTA=1.6*DELTA_original;
                if strcmp(models{modeli},'RBC_EZ') || strcmp(models{modeli},'RBC_EZ_adjCost')
                    DELTA=DELTA_original;
                elseif strcmp(models{modeli},'RBC_EZ_adjCost_Calvo')
                    DELTA=1.3*DELTA_original;
                end
            end
            addpath([main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order) '\files'])
            load([main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order) '\model'],'model')

            %%%%%%%
            n_x=model.n_x;
            n_x1=model.n_x1;
            n_x2=model.n_x2;
            n_y=model.n_y; 
            n_f=model.n_f;

            % build a Smolyak grid

            d=n_x;mu=order;
            Smol_grid = Smolyak_Grid(d,mu,model.Smolyak_elem_iso); 

            grid=Smol_grid';

            grid_points=size(grid,2);

            x_grid=repmat(cheb_c0,1,grid_points)+repmat(DELTA,1,grid_points).*grid; %unscaled grid

            % Compute the Smolyak coefficients that interpolate the initial
            % guess
            gh1=zeros(n_y+n_x1,grid_points);
            pert_coeffs=reshape(pert_coeffs,n_y+n_x1,[]); % use 3rd order perturbation as the initial guess
            [U2,W2]=create_UW(n_x,2);
            [U3,W3]=create_UW(n_x,3);
            temp=zeros(size(pert_coeffs,2),1);
            for i=1:grid_points
                kron1=x_grid(:,i)-nxss;
                temp(1:1+n_x)=[1;kron1];
                kron2=kron(kron1,kron1);
                temp(n_x+2:n_x+1+nchoosek(n_x+1,2))=W2*kron2;
                kron3=kron(kron2,kron1);
                temp(n_x+2+nchoosek(n_x+1,2):end)=W3*kron3;
                gh1(:,i)=pert_coeffs*temp;
            end

            Tcheb=Tcheb_fun(grid,[]);

            Tcoeffs=Tcheb'/gh1';
            Tcoeffs=Tcoeffs';

            coeffs=Tcoeffs(:); % this is the initial guess for the Smolyak coefficients

            model.maxload=smol_maxload;
            
            %  run first time to compile and update model indices
            [~,~,model]=loop_cheb(coeffs,model,params,x_grid(:,1),cheb_c0,DELTA,nep,P,eta,grid_maxload);

            % Solve with Newton method
            mkdir([main_folder '\' smol_folder '\' models{modeli} '\order' num2str(order)])

            try
                maxload=grid_maxload;
                fprintf(['\n' models{modeli} ': solving Smolyak, order ' num2str(order) '\n'])
                start=tic;
                [ ncoeffs,R,J,iter] = smol_newton( coeffs,model,params,x_grid,cheb_c0,DELTA,nep,P,eta,maxload,tolX,tolF,maxiter );
                time=toc(start);
                coeffs=ncoeffs;
            catch % newton failed
                coeffs=NaN*zeros(size(coeffs));
                R=[];
                J=[];
                iter=NaN;
                time=NaN;
            end
            save([main_folder '\' smol_folder '\' models{modeli} '\order' num2str(order)  '\solution_smol'],'coeffs','cheb_c0','DELTA','R','J','iter','time');
            rmpath([main_folder '\files_for_smolyak\' models{modeli} '\order' num2str(order) '\files'])
        end
    end
end

cd(main_folder)