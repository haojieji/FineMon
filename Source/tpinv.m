function pinvX = tpinv(X)

% tinv(X) is the inverse of the tensor X of size n*n*n3.
%   A warning message is printed if X is badly scaled or
%   nearly singular.
%
% version 1.0 - 14/06/2018
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

n = size(X);
n1 = n(1);
n2=n(2);
if length(n)==3
    n3=n(3);
else
    n3=1;
end

X = fft(X,[],3);
pinvX = zeros(n2,n1,n3);
for i=1:n3
    pinvX(:,:,i)=pinv(X(:,:,i));
end
pinvX = ifft(pinvX,[],3);