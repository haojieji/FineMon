%----- subroutine_ESTC : recover the un-sampled data by ESTC
%----- main inputs: Pu = H(U)*(H(U)^T_{\Omega}*H(U)_{\Omega}-\lambda*I)^-*H(u)_{\Omega},  Mi_omega 
%----- main outputs: the recovered data R
function [R,com]=subroutine_ESTC(R,Pu,Mi_omega,i,omega_slice_index,com)
    Mi_omega_1=Mi_omega(omega_slice_index,1,:);
    R(:,i,:)=tprod(Pu,Mi_omega_1);
    R(omega_slice_index,i,:)=Mi_omega(omega_slice_index,1,:);
    com(:,i)=1;
end

function R = TC(U_W_t,omega_slice_index,Mi_omega)
U_W_omega=U_W_t(omega_slice_index,:,:)
pesudoinverse_U_W_omega_t=tpinv(U_W_omega)
Pu=tprod(U_W_t,pesudoinverse_U_W_omega_t)

Mi_omega_1=Mi_omega(omega_slice_index,1,:);
R(:,i,:)=tprod(Pu,Mi_omega_1);
R(omega_slice_index,i,:)=Mi_omega(omega_slice_index,1,:);
end