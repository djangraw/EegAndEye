function [lambda_iter,sigmasq_iter] = GetLambdaIteratively(R,iEvents,new_filename)

% Created 8/27/14 by DJ.

MAX_ITER = 100;
CONVERGENCE_THRESHOLD = 1e-6;

% Update results struct
R = UpdateGlmResultsFormat(R);
UnpackStruct(R);
% Crop regressor events
regressor_events{iLevel} = regressor_events{iLevel}(iEvents);
% Overwrite filename
filenames{iLevel} = new_filename;

if ~exist('Y','var')
    fprintf('%s Getting matrices...',datestr(now,16))
    tic;    
    [X,Y] = GetGlmMatrices(EEG,regressor_events{iLevel},influence{iLevel},artifact_events,artifact_influence{iLevel},demeanX,normalizeX);
    t=toc;
    fprintf('took %.1f seconds.\n',t);
end   
% Normalize Matrices
X = double(X);
for i=1:size(X,2)
    X(:,i) = X(:,i)/sqrt(X(:,i)'*X(:,i));
end
for i=1:size(Y,2)
    Y(:,i) = Y(:,i)/sqrt(Y(:,i)'*Y(:,i));
end 

% Multiply matrices
if ~exist('XTY','var')
    fprintf('%s Multiplying matrices...',datestr(now,16))
    tic;
    XTX = X'*X;
    XTY = X'*Y;
    t=toc;
    fprintf('took %.1f seconds.\n',t);
end
% Solve
fprintf('%s Solving for betas...',datestr(now,16))
tic;
D = size(XTY,2);
p = size(XTX,2);

% Set up
lambda_iter = nan(1,MAX_ITER);
sigmasq_iter = nan(1,MAX_ITER);
beta = (XTX^(-1)*XTY)';
for j=1:MAX_ITER
    fprintf('Iter %d...',j);
    [sigmasq_all, sigmasq_elec] = GetGlmMse(X,Y,beta);
    k_elec = nan(1,D);
    for i=1:D
        k_elec(i) = p*sigmasq_elec(i)/sum(beta(i,:).^2);
    end
    lambda_iter(j) = mean(k_elec);
    sigmasq_iter(j) = sigmasq_all;    
    fprintf('o^2 = %g, lambda = %g\n',sigmasq_iter(j),lambda_iter(j));
    if j>1 && abs(lambda_iter(j)-lambda_iter(j-1))<CONVERGENCE_THRESHOLD
        disp('Converged! ')
        break
    elseif j==MAX_ITER
        disp('Failed to converge!')
    end
    beta = ( (XTX+lambda_iter(j)*eye(p))^(-1)*XTY )'; % size Dxp    
end
t=toc;
fprintf('took %.1f seconds.\n',t);

lambda_iter = lambda_iter(~isnan(lambda_iter));
sigmasq_iter = sigmasq_iter(~isnan(sigmasq_iter));

% Place converged lambda
lambda = {lambda_iter(end)};
% Place Response functions
responseFns{iLevel} = reshape(beta,[D,Nh,Nr]);       

% Save results
if ~isempty(filenames{iLevel})
    % Update status
    fprintf('%s Saving Results as %s...',datestr(now,16), filenames{iLevel});
    tic;
    % Save results
    save(filenames{iLevel},'X','Y','responseFns','tResponse','regressor_events',...
        'filenames','iLevel','artifact_events','artifact_influence','dataset','offset','influence','stddev',...
        'vthresh','method','trial_rej_rules','lambda','demeanX','normalizeX');
    t=toc;
    fprintf('took %.1f seconds.\n',t);
end
