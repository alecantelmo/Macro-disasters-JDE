% Discretize Gaussian shocks By Monomials

%%%%%%%%%% Monomial rules are taken from Judd, Maliar, Maliar and Valero 2014
n_e=size(eta_mat,2);
%%% 1st order monomial
[n_nodes,epsi_nodes,weight_nodes] = Monomials_1(n_e-1,eye(n_e-1));

nep1=epsi_nodes';
P1=weight_nodes';

% add the disaster shock
nep1=[zeros(1,size(nep1,2))-MUD,ones(1,size(nep1,2))-MUD;repmat(nep1,1,2)];
P1=[(1-prob_disaster)*P1,prob_disaster*P1];

%%% 2nd order monomial
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(n_e-1,eye(n_e-1));

nep2=epsi_nodes';
P2=weight_nodes';

% add the disaster shock
nep2=[zeros(1,size(nep2,2))-MUD,ones(1,size(nep2,2))-MUD;repmat(nep2,1,2)];
P2=[(1-prob_disaster)*P2,prob_disaster*P2];