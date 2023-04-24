%
% © Copyright, Oren Levintal, June 13, 2016.

if approx>=1
    nPhix=Phi_d_c1(nx(:)',params,model.Phi_indc);
    nPhix.tsize=[n_x2,repmat(n_x,1,1)]; % change size to account for the perturbation parameter
    nPhix=sptensor2spmat(unfold(nPhix));
end
if approx>=2
    nPhixxc=Phi_d_c2(nx(:)',params,model.Phi_indc);
    [nPhixx,model.pertind{pertindi}]=uncompressderivs( nPhixxc,2,n_x_tp,ones(n_x2,n_x_tp),model.pertind{pertindi} );
    pertindi=pertindi+1;
    nPhixx.tsize=[n_x2,repmat(n_x,1,2)]; % change size to account for the perturbation parameter
    nPhixx=sptensor2spmat(unfold(nPhixx));
end
if approx>=3
    nPhixxxc=Phi_d_c3(nx(:)',params,model.Phi_indc);
    [nPhixxx,model.pertind{pertindi}]=uncompressderivs( nPhixxxc,3,n_x_tp,ones(n_x2,n_x_tp),model.pertind{pertindi} );
    pertindi=pertindi+1;
    nPhixxx.tsize=[n_x2,repmat(n_x,1,3)]; % change size to account for the perturbation parameter
    nPhixxx=sptensor2spmat(unfold(nPhixxx));
end
if approx>=4
    nPhixxxxc=Phi_d_c4(nx(:)',params,model.Phi_indc);
    [nPhixxxx,model.pertind{pertindi}]=uncompressderivs( nPhixxxxc,4,n_x_tp,ones(n_x2,n_x_tp),model.pertind{pertindi} );
    pertindi=pertindi+1;
    nPhixxxx.tsize=[n_x2,repmat(n_x,1,4)]; % change size to account for the perturbation parameter
    nPhixxxx=sptensor2spmat(unfold(nPhixxxx));
end