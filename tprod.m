function C = tprod(A,B)

% Tensor-tensor product of two 3 way tensors: C = A*B
% A - n1*n2*n3 tensor
% B - n2*l*n3  tensor
% C - n1*l*n3  tensor
%
% version 2.0 - 09/10/2017
%
% Written by Canyi Lu (canyilu@gmail.com)
%
%
% References: 
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin and Shuicheng
% Yan, Tensor Robust Principal Component Analysis with A New Tensor Nuclear
% Norm, arXiv preprint arXiv:1804.03728, 2018
%

%[n1,n2,n3] = size(A);
n = size(A);
n1 = n(1);
n2=n(2);
if length(n)==3
    n3=n(3);
else
    n3=1;
end
%[m1,m2,m3] = size(B);
m = size(B);
m1 = m(1);
m2=m(2);
if length(m)==3
    m3=m(3);
else
    m3=1;
end

if n2 ~= m1 || n3 ~= m3 
    error('Inner tensor dimensions must agree.');
end

A = fft(double(A),[],3);
B = fft(double(B),[],3);
C = zeros(n1,m2,n3);

% first frontal slice
C(:,:,1) = A(:,:,1)*B(:,:,1);
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    C(:,:,i) = A(:,:,i)*B(:,:,i);
    C(:,:,n3+2-i) = conj(C(:,:,i));
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    C(:,:,i) = A(:,:,i)*B(:,:,i);
end
C = ifft(C,[],3);