function R_inv2=Find_Inverse_Matrix(R_inv,phi_f,psi_f,beta_f)
R_inv2=zeros(size(R_inv),length(phi_f));
for i=1:length(phi_f)
    R_inv2(:,:,i)=vpa(subs(R_inv));