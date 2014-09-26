n = 1000;
p = 100;
D = 10;


X = rand(n,p);
X(:,end) = X(:,1)+X(:,2);

beta = rand(p,D);

Y = X*beta;

%% Solve with OLS
beta0 = (X'*X)^-1*X'*Y;

subplot(2,2,1);
imagesc(beta-beta);
set(gca,'clim',[-1 1]);
colorbar
subplot(2,2,2);
imagesc(beta0-beta);
set(gca,'clim',[-1 1]);
colorbar

%% Use SVD
[U,S,V] = svds(X,p-1);
X1 = U*S;
beta1_pc = (X1'*X1)^-1*X1'*Y;
beta1 = V*beta1_pc;
subplot(2,2,3);
imagesc(beta1-beta);
set(gca,'clim',[-1 1]);
colorbar

%% Use PCA
[U1,S1,V1] = svd(X);
X1 = U1(:,1:p-1) * S1(1:p-1,1:p-1);
beta1_pc = (X1'*X1)^-1*X1'*Y;
beta1 = V1(:,1:p-1) * beta1_pc;
subplot(2,2,4);
imagesc(beta1-beta);
set(gca,'clim',[-1 1]);
colorbar

