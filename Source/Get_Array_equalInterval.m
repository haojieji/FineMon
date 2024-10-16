%{
equal interval sampling
input£ºthe number of rows n1, the number of samples m
output£ºarr,(1 if the row is sampled, 0 else)
%}
function arr=Get_Array_equalInterval(n1,m)
    interval = ceil(n1/m);
    Omega = 1:interval:n1;
    Omega_0=zeros(1,n1);
    Omega_0(Omega)=1;
    arr=Omega_0';
end