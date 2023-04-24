%%%%%%%%%%%%% code from tp to compute fv,fvv,fvvv,fvvvv, where v is the
%%%%%%%%%%%%% vector of variables of the tp system (the pert system has a
%%%%%%%%%%%%% larger v).
%
% © Copyright, Oren Levintal, June 13, 2016.

nv=[zeros(n_y,1);nyss;nxss;nxss]; % variables of tp

n_u=model.n_u;
nu=zeros(n_u,1);
npreu=preu_fun(nv(model.preuvars),params); % all predetermined u
nu(model.preurows)=npreu;
nz=[nv;nu]; % auxiliary variables
pref=pretilf_fun(nz(model.pretilfzvars),params); % since i already have z, i use pretilf to evaluate f.


% compute derivatives of prePI and pretilf w.r.t z (prePI=nonstochastic rows of PI, pretilf=nonstochastic rows of f)

prePIz_full=prePI_d1(nz',params,model.prePI_ind_u);
prefvars=model.prefvars;
n_prefvars=length(prefvars);
pretilfz=pretilf_tilf_d1([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
if approx>=2
    prePIzz_full=prePI_d2(nz',params,model.prePI_ind_u);
    pretilfzz=pretilf_tilf_d2([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
end
if approx>=3
    prePIzzz_full=prePI_d3(nz',params,model.prePI_ind_u);
    pretilfzzz=pretilf_tilf_d3([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
end
if approx>=4
    prePIzzzz_full=prePI_d4(nz',params,model.prePI_ind_u);
    pretilfzzzz=pretilf_tilf_d4([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
end

% extract derivatives of prePI w.r.t v from prePIz, prePIzz, ...
if ~isfield(model,'ind')
    model.ind=cell(model.totindi,1);
end
indi=5;

if approx>=1
    [prePIz,model.ind{indi}]=extract(prePIz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIv,model.ind{indi}]=extract(prePIz_full,model.prefzvars,model.prefvars,0,model.ind{indi});
    indi=indi+1;
end
if approx>=2
    [prePIzz,model.ind{indi}]=extract(prePIzz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIvvc,model.ind{indi}]=extract(prePIzz_full,model.prefzvars,model.prefvars,1,model.ind{indi});
    indi=indi+1;
end
if approx>=3
    [prePIzzz,model.ind{indi}]=extract(prePIzzz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIvvvc,model.ind{indi}]=extract(prePIzzz_full,model.prefzvars,model.prefvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'prefvars_chain3c_M2')
        [ tempM ] = chainsM( n_prefvars,3 );
        model.prefvars_chain3c_M2=tempM{2};
        clear tempM
    end
end
if approx>=4
    [prePIzzzz,model.ind{indi}]=extract(prePIzzzz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIvvvvc,model.ind{indi}]=extract(prePIzzzz_full,model.prefzvars,model.prefvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'prefvars_chain4c_M2')
        [ tempM ] = chainsM( n_prefvars,4 );
        model.prefvars_chain4c_M2=tempM{2};
        model.prefvars_chain4c_M3=tempM{3};
        model.prefvars_chain4c_M4=tempM{4};
        clear tempM
    end
end


% Use high order chain rules to compute derivatives of prePI w.r.t v
totindi=5+8;
if approx==1
    indi=totindi;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
elseif approx==2
    indi=totindi+(model.pre_n-1)*(1)+1;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=chain2c_tensor(prePIz,prePIzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
    indi=indi+1;

    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvc,model.ind{indi}]=chain2c_tensor(pretilfz,pretilfzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [ prefvv,model.ind{indi} ] = uncompressderivs( prefvvc,2,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

elseif approx==3
    indi=totindi+(model.pre_n-1)*(1+3)+1+4;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
        indi=indi+1;
        [prePIvvvc,model.ind{indi}]=chain3c_tensor(prePIz,prePIzz,prePIzzz,...
            prePIv,prePIvvc,prePIvvvc,model.prefvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=chain2c_tensor(prePIz,prePIzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
    indi=indi+1;

    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvc,model.ind{indi}]=chain2c_tensor(pretilfz,pretilfzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvvc,model.ind{indi}]=chain3c_tensor(pretilfz,pretilfzz,pretilfzzz,...
        prePIv,prePIvvc,prePIvvvc,...
        model.prefvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [ prefvv,model.ind{indi} ] = uncompressderivs( prefvvc,2,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

    [ prefvvv,model.ind{indi} ] = uncompressderivs( prefvvvc,3,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

elseif approx==4
    indi=totindi+(model.pre_n-1)*(1+3+4)+1+4+6;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=colsort(prePIvvc,model.ind{indi});
        indi=indi+1;
        [prePIvvvvc,model.ind{indi}]=chain4c_tensor(prePIz,prePIzz,prePIzzz,prePIzzzz,...
            prePIv,prePIvvc,prePIvvvc,prePIvvvvc,...
            model.prefvars_chain4c_M2,model.prefvars_chain4c_M3,model.prefvars_chain4c_M4,...
            model.ind{indi},n_ind,maxload,'vec',model.prezz,model.maxzz);
        indi=indi+1;
        [prePIvvvc,model.ind{indi}]=chain3c_tensor(prePIz,prePIzz,prePIzzz,...
            prePIv,prePIvvc,prePIvvvc,model.prefvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=chain2c_tensor(prePIz,prePIzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
    indi=indi+1;
    [prePIvvc,model.ind{indi}]=colsort(prePIvvc,model.ind{indi});
    indi=indi+1;

    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [prefvvc,model.ind{indi}]=chain2c_tensor(pretilfz,pretilfzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvvc,model.ind{indi}]=chain3c_tensor(pretilfz,pretilfzz,pretilfzzz,...
        prePIv,prePIvvc,prePIvvvc,...
        model.prefvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvvvc,model.ind{indi}]=chain4c_tensor(pretilfz,pretilfzz,pretilfzzz,pretilfzzzz,...
        prePIv,prePIvvc,prePIvvvc,prePIvvvvc,...
        model.prefvars_chain4c_M2,model.prefvars_chain4c_M3,model.prefvars_chain4c_M4,...
        model.ind{indi},n_ind,maxload,'vec',model.pretilfz,model.maxtilfz);
    indi=indi+1;

    [ prefvv,model.ind{indi} ] = uncompressderivs( prefvvc,2,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

    [ prefvvv,model.ind{indi} ] = uncompressderivs( prefvvvc,3,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

    [ prefvvvv,model.ind{indi} ] = uncompressderivs( prefvvvvc,4,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;
end
totindi=totindi+(model.pre_n-1)*(1+3+4+6)+1+4+6+9;

%%%%%%%%%%%%% now stochastic equations
n_nodes=1;
nv_vec=[nyss;nyss;nxss;nxss];
stochfvars=model.stochfvars;
n_stochfvars=length(stochfvars);
nstochv_vec=nv_vec(stochfvars,:); % stochastic vars.

nstochu_vec=stochu_fun(nv_vec(model.stochuvars,:),params); % all stochastic u vars
nu_vec=repmat(nu,1,n_nodes);
nu_vec(model.stochurows,:)=nstochu_vec;
nz_vec=[nv_vec;nu_vec];
stochf_vec=stochtilf_fun(nz_vec(model.stochtilfzvars,:),params);

% compute derivatives of stochPI and stochtilf w.r.t z

stochPIz_vec=stochPI_d1(nz_vec',params,model.stochPI_ind_u);
prePIz_full.vals=repmat(prePIz_full.vals,n_nodes,1);

i1=tfind(prePIz_full);
i2=tfind(stochPIz_vec);
PIz_vec=sptensor([i1;i2],[prePIz_full.cols;stochPIz_vec.cols],[prePIz_full.vals,stochPIz_vec.vals],n_z,[n_z]);
stochtilfz_vec=stochtilf_tilf_d1([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);

if approx>=2
    stochPIzz_vec=stochPI_d2(nz_vec',params,model.stochPI_ind_u);
    prePIzz_full.vals=repmat(prePIzz_full.vals,n_nodes,1);

    i1=tfind(prePIzz_full);
    i2=tfind(stochPIzz_vec);
    PIzz_vec=sptensor([i1;i2],[prePIzz_full.cols;stochPIzz_vec.cols],[prePIzz_full.vals,stochPIzz_vec.vals],n_z,[n_z,n_z]);
    stochtilfzz_vec=stochtilf_tilf_d2([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);
end
if approx>=3
    stochPIzzz_vec=stochPI_d3(nz_vec',params,model.stochPI_ind_u);
    prePIzzz_full.vals=repmat(prePIzzz_full.vals,n_nodes,1);

    i1=tfind(prePIzzz_full);
    i2=tfind(stochPIzzz_vec);
    PIzzz_vec=sptensor([i1;i2],[prePIzzz_full.cols;stochPIzzz_vec.cols],[prePIzzz_full.vals,stochPIzzz_vec.vals],n_z,[n_z,n_z,n_z]);
    stochtilfzzz_vec=stochtilf_tilf_d3([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);
end
if approx>=4
    stochPIzzzz_vec=stochPI_d4(nz_vec',params,model.stochPI_ind_u);
    prePIzzzz_full.vals=repmat(prePIzzzz_full.vals,n_nodes,1);

    i1=tfind(prePIzzzz_full);
    i2=tfind(stochPIzzzz_vec);
    PIzzzz_vec=sptensor([i1;i2],[prePIzzzz_full.cols;stochPIzzzz_vec.cols],[prePIzzzz_full.vals,stochPIzzzz_vec.vals],n_z,[n_z,n_z,n_z,n_z]);
    stochtilfzzzz_vec=stochtilf_tilf_d4([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);
end

% extract PIv,Pivv,... from PIz,PIzz,...

totindi=totindi+3;
totindi=totindi+4;
totindi=totindi+2;
totindi=totindi+1;
totindi=totindi+1;
totindi=totindi+9;


indi=totindi;
if approx>=1
    [stochPIz_vec,model.ind{indi}]=extract(PIz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIv_vec,model.ind{indi}]=extract(PIz_vec,model.stochfzvars,model.stochfvars,0,model.ind{indi});
    indi=indi+1;  
end
if approx>=2
    [stochPIzz_vec,model.ind{indi}]=extract(PIzz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIvvc_vec,model.ind{indi}]=extract(PIzz_vec,model.stochfzvars,model.stochfvars,1,model.ind{indi});
    indi=indi+1;
end
if approx>=3
    [stochPIzzz_vec,model.ind{indi}]=extract(PIzzz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIvvvc_vec,model.ind{indi}]=extract(PIzzz_vec,model.stochfzvars,model.stochfvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'stochfvars_chain3c_M2')
        [ tempM ] = chainsM( n_stochfvars,3 );
        model.stochfvars_chain3c_M2=tempM{2};
        clear tempM
    end
end
if approx>=4
    [stochPIzzzz_vec,model.ind{indi}]=extract(PIzzzz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIvvvvc_vec,model.ind{indi}]=extract(PIzzzz_vec,model.stochfzvars,model.stochfvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'stochfvars_chain4c_M2')
        [ tempM ] = chainsM( n_stochfvars,4 );
        model.stochfvars_chain4c_M2=tempM{2};
        model.stochfvars_chain4c_M3=tempM{3};
        model.stochfvars_chain4c_M4=tempM{4};
        clear tempM
    end
end

totindi=totindi+8;

% compute derivatives of stochf w.r.t v by high order chain rules

if approx==1
    indi=totindi;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

elseif approx==2
    indi=totindi+(model.stoch_n-1)*(1)+1;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=chain2c_tensor(stochPIz_vec,stochPIzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
    indi=indi+1;

    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvc_vec,model.ind{indi}]=chain2c_tensor(stochtilfz_vec,stochtilfzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [ stochfvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvc_vec,2,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvc_vec
    
elseif approx==3
    indi=totindi+(model.stoch_n-1)*(1+3)+1+4;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvvc_vec,model.ind{indi}]=chain3c_tensor(stochPIz_vec,stochPIzz_vec,stochPIzzz_vec,...
            stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,model.stochfvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=chain2c_tensor(stochPIz_vec,stochPIzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
    indi=indi+1;

    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvc_vec,model.ind{indi}]=chain2c_tensor(stochtilfz_vec,stochtilfzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvvc_vec,model.ind{indi}]=chain3c_tensor(stochtilfz_vec,stochtilfzz_vec,stochtilfzzz_vec,...
        stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,...
        model.stochfvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [ stochfvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvc_vec,2,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvc_vec

    [ stochfvvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvvc_vec,3,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvvc_vec

elseif approx==4
    indi=totindi+(model.stoch_n-1)*(1+3+4)+1+4+6;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=colsort(stochPIvvc_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvvvc_vec,model.ind{indi}]=chain4c_tensor(stochPIz_vec,stochPIzz_vec,stochPIzzz_vec,stochPIzzzz_vec,...
            stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,stochPIvvvvc_vec,...
            model.stochfvars_chain4c_M2,model.stochfvars_chain4c_M3,model.stochfvars_chain4c_M4,...
            model.ind{indi},n_ind,maxload,'vec',model.stochzz,model.maxzz);
        indi=indi+1;
        [stochPIvvvc_vec,model.ind{indi}]=chain3c_tensor(stochPIz_vec,stochPIzz_vec,stochPIzzz_vec,...
            stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,model.stochfvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=chain2c_tensor(stochPIz_vec,stochPIzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
    indi=indi+1;
    [stochPIvvc_vec,model.ind{indi}]=colsort(stochPIvvc_vec,model.ind{indi});
    indi=indi+1;

    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [stochfvvc_vec,model.ind{indi}]=chain2c_tensor(stochtilfz_vec,stochtilfzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvvc_vec,model.ind{indi}]=chain3c_tensor(stochtilfz_vec,stochtilfzz_vec,stochtilfzzz_vec,...
        stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,...
        model.stochfvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvvvc_vec,model.ind{indi}]=chain4c_tensor(stochtilfz_vec,stochtilfzz_vec,stochtilfzzz_vec,stochtilfzzzz_vec,...
        stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,stochPIvvvvc_vec,...
        model.stochfvars_chain4c_M2,model.stochfvars_chain4c_M3,model.stochfvars_chain4c_M4,...
        model.ind{indi},n_ind,maxload,'vec',model.stochtilfz,model.maxtilfz);
    indi=indi+1;

    [ stochfvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvc_vec,2,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvc_vec

    [ stochfvvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvvc_vec,3,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvvc_vec

    [ stochfvvvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvvvc_vec,4,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvvvc_vec

end

%%%%%%%%%%%%% convert to sparse matrices

if approx>=1
    prefv=changecols(prefv,model.prefvars,n_v_tp,1);
    tempv=1:n_v;
    tempv(2*n_y+n_x)=[];
    tempv(end)=[];
    prefv=changecols(prefv,tempv,n_v,1);
    prefv=changerows(prefv,model.prefrows,n_f);
    
    stochfv_vec=changecols(stochfv_vec,model.stochfvars,n_v_tp,1);
    stochfv_vec=changecols(stochfv_vec,tempv,n_v,1);
    stochfv_vec=changerows(stochfv_vec,model.stochfrows,n_f);

    [prei,prej,prevals]=tfind(unfold(prefv));
    [stochi,stochj,stochvals]=tfind(unfold(stochfv_vec));
    
    fv=sparse(double([prei;stochi]),double([prej;stochj]),double([prevals(:);stochvals(:)]),n_f,n_v);
end
if approx>=2
    prefvv=changecols(prefvv,model.prefvars,n_v_tp,1);
    prefvv=changecols(prefvv,model.prefvars,n_v_tp,2);
    prefvv=changecols(prefvv,tempv,n_v,1);
    prefvv=changecols(prefvv,tempv,n_v,2);
    prefvv=changerows(prefvv,model.prefrows,n_f);
    
    stochfvv_vec=changecols(stochfvv_vec,model.stochfvars,n_v_tp,1);
    stochfvv_vec=changecols(stochfvv_vec,model.stochfvars,n_v_tp,2);
    stochfvv_vec=changecols(stochfvv_vec,tempv,n_v,1);
    stochfvv_vec=changecols(stochfvv_vec,tempv,n_v,2);
    stochfvv_vec=changerows(stochfvv_vec,model.stochfrows,n_f);

    [prei,prej,prevals]=tfind(unfold(prefvv));
    [stochi,stochj,stochvals]=tfind(unfold(stochfvv_vec));
    
    alli=double([prei;stochi]);
    allj=double([prej;stochj]);
    allvals=double([prevals(:);stochvals(:)]);
    
    newi=sub2ind([n_f,n_v^2],alli,allj);
    newj=ones(size(newi));
    
    fvv=sparse(newi,newj,allvals,n_f*n_v^2,1);
end
if approx>=3
    prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,1);
    prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,2);
    prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,3);
    prefvvv=changecols(prefvvv,tempv,n_v,1);
    prefvvv=changecols(prefvvv,tempv,n_v,2);
    prefvvv=changecols(prefvvv,tempv,n_v,3);
    prefvvv=changerows(prefvvv,model.prefrows,n_f);
    
    stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,1);
    stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,2);
    stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,3);
    stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,1);
    stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,2);
    stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,3);
    stochfvvv_vec=changerows(stochfvvv_vec,model.stochfrows,n_f);
    
    [prei,prej,prevals]=tfind(unfold(prefvvv));
    [stochi,stochj,stochvals]=tfind(unfold(stochfvvv_vec));
    
    alli=double([prei;stochi]);
    allj=double([prej;stochj]);
    allvals=double([prevals(:);stochvals(:)]);
    
    newi=sub2ind([n_f,n_v^3],alli,allj);
    newj=ones(size(newi));
    
    fvvv=sparse(newi,newj,allvals,n_f*n_v^3,1);
end
if approx>=4
    prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,1);
    prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,2);
    prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,3);
    prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,4);
    prefvvvv=changecols(prefvvvv,tempv,n_v,1);
    prefvvvv=changecols(prefvvvv,tempv,n_v,2);
    prefvvvv=changecols(prefvvvv,tempv,n_v,3);
    prefvvvv=changecols(prefvvvv,tempv,n_v,4);
    prefvvvv=changerows(prefvvvv,model.prefrows,n_f);
    
    stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,1);
    stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,2);
    stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,3);
    stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,4);
    stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,1);
    stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,2);
    stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,3);
    stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,4);
    stochfvvvv_vec=changerows(stochfvvvv_vec,model.stochfrows,n_f);
    
    [prei,prej,prevals]=tfind(unfold(prefvvvv));
    [stochi,stochj,stochvals]=tfind(unfold(stochfvvvv_vec));
    
    alli=double([prei;stochi]);
    allj=double([prej;stochj]);
    allvals=double([prevals(:);stochvals(:)]);
    
    newi=sub2ind([n_f,n_v^4],alli,allj);
    newj=ones(size(newi));
    
    fvvvv=sparse(newi,newj,allvals,n_f*n_v^4,1);
end