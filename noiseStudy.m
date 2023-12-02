% noise study

clear; clc;
r=4;
A = tprod(rand(100,r,100), rand(r,100,100));
[U,S,V] = tsvd(A);
% 计算秩 r=非零奇异tube
r1=0;
for i = 1:size(S,1)
    if find(S(i,i,:)>1e-10, 1)
        r1=r1+1;
    end
end

%计算奇异tube的范数
engery = zeros(1,r1);
for i = 1:r1
    tubei(:) = S(i,i,:);
    engery(i) = norm(tubei);
end

% 高斯噪声E
% B = tprod(rand(10,2,10), rand(2,10,10));
% E = mean(A) + engery(1)*B;
% E = B;

% amount
p = 0.05;
% perm = randperm(size(A,1)*size(A,2)*size(A,3));
% samplenumber = p*length(perm);
% Omega = zeros(size(A));
% Omega(perm(1:samplenumber)) = 1;

perm = randperm(size(A,1)*size(A,2));
samplenumber = p*length(perm);

% Omega = zeros(size(A));
% for i = 1:size(Omega,3)
%     Omegai = zeros(size(Omega,1),size(Omega,2));
%     Omegai(perm(1:samplenumber)) = 1;
%     Omega(:,:,i) = Omegai;
% end

% A_E = A + E.*Omega;
A_E = A;
% A_E(perm(1:samplenumber)) = A_E(perm(1:samplenumber)) +  rand(1,samplenumber);

omega = zeros(size(A,1),size(A,2));
omega(perm(1:samplenumber)) = 1;
for i=1:size(omega,1)
    for j=1:size(omega,2)
        if omega(i,j) == 1
            A_E(i,j,:) =  A_E(i,j,:) + rand(1);
        end
    end
end

[U_E,S_E,V_E] = tsvd(A_E);
% 计算秩 r=非零奇异tube
r2=0;
for i = 1:size(S_E,1)
    if find(S_E(i,i,:)>1e-10, 1)
        r2=r2+1;
    end
end

%%
clear
Ps = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
results_r = zeros(1,length(Ps));
r=10;
A = tprod(rand(100,r,100), rand(r,100,100));
[U,S,V] = tsvd(A);
% 计算秩 r=非零奇异tube
r1=0;
for i = 1:size(S,1)
    if find(S(i,i,:)>1e-10, 1)
        r1=r1+1;
    end
end

%计算奇异tube的范数
engery = zeros(1,r1);
for i = 1:r1
    tubei(:) = S(i,i,:);
    engery(i) = norm(tubei);
end

for i =  1:length(Ps)
    p = Ps(i);
    perm = randperm(size(A,1)*size(A,2));
    samplenumber = p*size(A,1)*size(A,2);
    A_E = A;
    omega = zeros(size(A,1),size(A,2));
    omega(perm(1:samplenumber)) = 1;
    for ii=1:size(omega,1)
        for jj=1:size(omega,2)
            if omega(ii,jj) == 1
                A_E(ii,jj,:) =  A_E(ii,jj,:) + 100*rand(1);
            end
        end
    end
    
    [U_E,S_E,V_E] = tsvd(A_E);
    % 计算秩 r=非零奇异tube
    r2=0;
    for j = 1:size(S_E,1)
        if find(S_E(j,j,:)>1e-10, 1)
            r2=r2+1;
        end
    end

    results_r(i) = r2;
end

figure
plot(Ps(1:10), results_r(1:10),  'o-','linewidth',5,'MarkerSize',5);
xticks(Ps(1:10))
yticks(10:20:100)
xlabel('Amount of noise');
ylabel('Interupted rank');
set(gca,'FontName','Times New Roman' ,'FontSize',24, 'FontWeight','bold');

%%
%Q: 一个线性相关列的噪声量为多少时, 该相关列转变为独立列
clear
r=2;
A = tprod(rand(100,r,100), rand(r,100,100));
[U,S,V] = tsvd(A);
% 计算秩 r=非零奇异tube
r1=getTrank(S);
rs=[];
rs=[rs r1];

idx = randperm(size(A,1));
for i = 1:size(A,1)
    %干扰第10列的i个位置
    A10 = A(:,1,:);
    E = zeros(i,1,size(A,3));
    % e = randn(i,1); %i个位置的噪声；对tube中的所有元素干扰相同的值
    for k = 1:size(A,3)
        e = randn(i,1); %对tube不同维度干扰不同值
        E(:,1,k) =  e;
    end
    A10(idx(1:i),:,:) = A10(idx(1:i),:,:) + E;
    A(:,1,:)=A10;
    %计算干扰后的秩
    [~,Si,~]=tsvd(A);
    r2=getTrank(Si);
    rs=[rs r2];
end



%%
%定理验证：对于一个低秩r的张量，噪音干扰的列越多，秩越高；当噪音干扰的slice个数超过(TM)^(1/4),秩满秩
clear;
r = 2;
A = tprod(rand(100,r,100), rand(r,100,100));
[U,S,V] = tsvd(A);
% 计算秩 r=非零奇异tube
r1=getTrank(S);
rs=[];

% 对列干扰
perm = randperm(size(A,2));
c=3;
bound = ceil(power(c*100*100, 1/4))-2;
for i = 1:bound
    AE=A;
    idx = perm(1:i);
    E = zeros(size(AE,1),length(idx),size(AE,3));
    for j = 1:i
        e = randn(size(A,1),1);
        for k = 1:size(A,3)
            E(:,j,k) = e;
        end
    end
    AE(:,idx,:) = AE(:,idx,:) + E;

    [~,SE,~] = tsvd(AE);
    r2 = getTrank(SE);

    rs = [rs r2];
end
samples = c*power(rs, 5/2).*log(2*rs);
samples_bound_isright = find(samples>100*100);
disp(samples_bound_isright);

%%
%定理验证：逐渐增加干扰的slice个数，看slice个数达到多少时，频率为1
% 对于一个低秩r的张量，噪音干扰的列越多，秩越高；当噪音干扰的slice个数超过(TM)^(1/4),秩满秩
clear;
r = 2;
A = tprod(rand(100,r,100), rand(r,100,100));
[U,S,V] = tsvd(A);
% 计算秩 r=非零奇异tube
r1=getTrank(S);
rs=[];
samples=[];
c = 24;
% 对列干扰
perm = randperm(size(A,2));
for i = 1:20
    % 每次干扰 i 个slice
    AE=A;
    idx = perm(1:i);
    E = zeros(size(AE,1),length(idx),size(AE,3));
    for j = 1:i
        e = randn(size(A,1),1);
        for k = 1:size(A,3)
            E(:,j,k) = e;
        end
    end
    AE(:,idx,:) = AE(:,idx,:) + E;
    [~,SE,~] = tsvd(AE);
    r2 = getTrank(SE);
    samples = [ samples power(r2, 5/2).*log(2*r2)];
    rs = [rs r2];
end

figure
plot(1:20, 6*samples,  'o-','linewidth',5,'MarkerSize',5);
hold on
plot(0:20, ones(1,21)*10000,  'k--','linewidth',4,'MarkerSize',5);
% xticks(Ps(1:10))
% yticks(0:2000:10000)
% ylim([0,10000])
xlabel('Amount of noise');
ylabel('Samples number');
set(gca,'FontName','Times New Roman' ,'FontSize',24, 'FontWeight','bold');

figure
plot(1:20, rs,  'o-','linewidth',5,'MarkerSize',5);
xlabel('Amount of noise');
ylabel('Interrupted rank');
set(gca,'FontName','Times New Roman' ,'FontSize',24, 'FontWeight','bold');

% 计算 tube 秩
function r1=getTrank(S)
    r1=0;
    for i = 1:size(S,1)
        if find(S(i,i,:)>1e-10, 1)
            r1=r1+1;
        end
    end
end