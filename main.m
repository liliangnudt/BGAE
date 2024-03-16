function [res] = main(X,Y,d,eta,lambda,beta,m)
rng_num = 0;
stream = RandStream.getGlobalStream;
reset(stream);

n = size(Y,1);
v = length(X);
c = length(unique(Y));

flag = 1;
iter = 0;
maxIter = 50;
IterMax = 50;
obj = [];

XX = [];
for p = 1 : v
    XX = [XX;X{p}];
    dim(p) = size(X{p},1);
end

for p = 1 : v
    Up{p} = zeros(dim(p), d);
    [IDXU,~]= kmeans(X{p}, d, 'MaxIter',200,'Replicates',30);
    for i = 1 : dim(p)
        U{i}(i,IDXU(i)) = 1;
    end
end

[XUV,~,~]=svds(XX', d);
V = XUV';

[XZ,~,~]=svds(XX', m);
[IDXZ,~] = kmeans(XZ, m, 'MaxIter',200,'Replicates',30);

for i = 1:n
    Z(i,IDXZ(i)) = 1;
end
Z = Z/(m) + (m-1)/m/m;

[XA,~,~]=svds(XX', d);
[~,A_temp] = kmeans(XA, m, 'MaxIter',200,'Replicates',30);
A = A_temp';

for p = 1 : v
    Wp{p} = eye(dim(p), m);
end

alpha = ones(1,v)/v;
gamma = ones(1,v)/v;

for p = 1 : v
    AAp{p} = zeros(size(X{p}));
end

while flag
    iter = iter + 1;
    %%
    for p = 1 : v
        H{p} = X{p} - Up{p} * V + (1 / beta) * AAp{p};
        Ep{p} = L21(H{p}, alpha(p)^2 / beta);
    end
    %%
    for p = 1 : v
        Up{p} = (X{p} - Ep{p} + 1 / beta * AAp{p}) * V';
    end
    %%
    d_Dcol = Z' * ones(n,1);
    ind = find(d_Dcol<1e-8);
    Z_2 = Z;
    Z_2(:,ind) = rand(n,length(ind));
    d_Dcol = Z_2' * ones(n,1);
    clear ind;
    A = V * Z_2 / diag(d_Dcol);
    %%
    for p = 1 : v
        WB = X{p} * Z;
        [UWB,~,VWB] = svd(WB,'econ');
        Wp{p} = UWB * VWB';
    end
    %%
    Z_temp = 0;
    for p = 1 : v
        Z_temp = Z_temp + gamma(p)^2 * X{p}' * Wp{p};
    end
    sum_gamma = 2 * lambda * norm(gamma,2)^2;
    et = -(eta * repmat((diag(A' * A))', n, 1) - 2 * (eta * V' * A + lambda * Z_temp));
    G = zeros(n,m); 
    G(1:m,:) = eye(m);
    [PreY,~,Z,~,~] = rank(et, G, c, sum_gamma, IterMax);
    %%
    Q_temp = 0;
    for p = 1 : v
        Qp{p} = X{p} - Ep{p} + 1 / beta * AAp{p};
        Q_temp = Q_temp + Qp{p}' * Up{p};
    end
    M = 2 * eta * Z * A' + beta * Q_temp;
    [UV,~,VV] = svd(M,'econ');
    V = (UV * VV')';
    %%
    Ma = zeros(v,1);
    for p = 1 : v
        al_temp1 = (X{p} - Up{p} * V)';
        al_temp2 = sqrt(sum(al_temp1.*al_temp1, 2));
        Ma(p) = sum(al_temp2);
    end
    Mafra = Ma.^(-1);
    Qa = 1/sum(Mafra);
    alpha = Qa*Mafra;
    %%
    Mb = zeros(v,1);
    for p = 1 : v
        Mb(p) = norm(Wp{p} * Z' - X{p},'fro')^2;
    end
    Mbfra = Mb.^(-1);
    Qb = 1/sum(Mbfra);
    gamma = Qb*Mbfra;
    %%
    for p = 1 : v
        AAp{p} = AAp{p} + beta * (X{p} - Up{p} * V - Ep{p});
    end
    beta = beta * 2;
    %%
    [obj(end+1)] = callobj(v, X, Up, V, A, Z, Wp, alpha, gamma, eta, lambda);
    if (iter>9) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>30 || obj(iter) < 1e-10)
        flag = 0;
    end
end
res = Clustering8Measure(Y, PreY);