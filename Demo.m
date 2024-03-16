clear;
clc;

addpath(genpath('./'));
Dataset_Path = './';
load(strcat(Dataset_Path,'Dermatology.mat'));

%%
n = size(Y,1);
v = length(X);
c = length(unique(Y));

for p = 1:v
    X{p} = mapstd(X{p}',0,1);
    X_dim(p) = size(X{p},1);
    d_min = min(X_dim);
end

m_set = [1]*c;
d_set = [1]*c;
eta = 1;
lambda = 1;
beta = 10;

%%
for AI = 1:length(m_set)
    m = m_set(AI);
    if m > d_min | m > n
        continue
    end
    for FI = 1:length(d_set)
        d = d_set(FI);
        if d > n | d > d_min
            continue
        end
        [res] = main(X,Y,d,eta,lambda,beta,m);
        fprintf('m:%2.0f \t d:%2.0f \t ACC:%4.2f \t NMI:%4.2f \t Pur:%4.2f \t Fscore:%4.2f \n',...
            [m d res(1)*100 res(2)*100 res(3)*100 res(4)*100]);
    end
end

