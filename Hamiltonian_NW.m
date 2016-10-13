function H = Hamiltonian_NW(Sites,layers,t,mu,alpha,V,delta)

sigma_x = [0,1;1,0]; sigma_y = [0,-1i;1i,0]; sigma_z = [1,0;0,-1];

h_v = V(1)*sigma_x + V(2)*sigma_y + V(3)*sigma_z;
h_onsite = [2*t - mu,0;0,2*t - mu] + h_v;
h_tx= [-t,alpha;-alpha,-t];
h_ty= [-t,-1i*alpha;-1i*alpha,-t];
% h_ty= [-0,0;0,0];
h_sc = [0,delta;-delta,0];

%====  Nambu representation

H00= [h_onsite,h_sc;h_sc',-conj(h_onsite)];
Htx = [h_tx,zeros(2);zeros(2),-conj(h_tx)];
Hty = [h_ty,zeros(2);zeros(2),-conj(h_ty)];

H0=kron(eye(layers),H00)+kron(diag(ones(1,layers-1),1),Hty)+kron(diag(ones(1,layers-1),-1),Hty');

Ht=kron(eye(layers),Htx);

H = kron(eye(Sites),H0)+kron(diag(ones(1,Sites-1),1),Ht)+kron(diag(ones(1,Sites-1),-1),Ht');

