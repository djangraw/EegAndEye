clear t_hat t_mvr rms_err_hat rms_err_mvr
for i=1:3
    fprintf('iter %d...\n',i);
    D = 50; % electrodes
    P = 8; % regressor types
    T = 75; % regressors
    N = 10000; % samples in experiment

    beta_true = rand(P*T,D);       
    X = zeros(N,P*T);
    for j=1:P
        iCols = (j-1)*T + (1:T);
        for k=1:100;
            iSpot = round((N-T-1)*rand(1));
            X(iSpot+(1:T),iCols) = X(iSpot+(1:T),iCols) + eye(T);
        end
    end
        
        
    Y = X*beta_true + 0.5*rand(N,D);

    %%
    tic
    X_inv = (X'*X)^(-1);
    beta_hat = X_inv*X'*Y;
    fprintf('least squares: ')
%     [beta_hat2,tH,Efit,Eres] = SimpleGLM(Y,[],X,[0 T-1],'leastsquares');
    disp('done!')
    
    t_hat(i) = toc;
    %%
    beta_mvr = zeros(size(beta_true));
    tic
    fprintf('mvregress: ')
    for j=1:D
        fprintf('chan %d/%d...',j,D);
        beta_mvr(:,j) = mvregress(X,Y(:,j));         
%         [beta_temp,tH,Efit,Eres] = SimpleGLM(Y(:,j),[],X,[0 T-1],'mvregress');
%         beta_mvr(:,j) = resize(beta_temp,P*T,1);
%         disp('done!')
    end    
    t_mvr(i) = toc;

    %% Get error
    rms_err_hat(i) = sqrt(mean((beta_true(:)-beta_hat(:)).^2));
    rms_err_mvr(i) = sqrt(mean((beta_true(:)-beta_mvr(:)).^2));
end

%% Print results
fprintf('Inverse Method: t=%g, err=%g\n',mean(t_hat),mean(rms_err_hat));
fprintf('MVRegress Method: t=%g, err=%g\n',mean(t_mvr),mean(rms_err_mvr));

%%
figure(101); clf;
subplot(3,3,1);
imagesc(beta_true); title('beta')
colorbar;
subplot(3,3,2);
imagesc(X); title('X');
subplot(3,3,3);
imagesc(Y); title('Y');

subplot(3,3,4);
imagesc(beta_hat); title('beta-hat');
colorbar;
subplot(3,3,5);
imagesc(beta_mvr); title('beta-mvregress');
colorbar;
subplot(3,3,6);
imagesc(beta_mvr-beta_hat); title('mvregress - hat');
colorbar;

subplot(3,3,7);
imagesc(beta_true-beta_hat); title('beta-hat error')
colorbar;
subplot(3,3,8);
imagesc(beta_true-beta_mvr); title('beta-mvregress error')
colorbar;