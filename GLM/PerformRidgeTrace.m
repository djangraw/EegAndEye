function PerformRidgeTrace(R,lambda_list,new_filename)

% PerformRidgeTrace(R,lambda_list)
%
% Created 3/19/14 by DJ.
% Updated 8/14/14 by DJ - added lambda input to RunEegGlm.

% Update results struct
R = UpdateGlmResultsFormat(R);

UnpackStruct(R);
nLevels = length(regressor_events);

% overwrite lambda & filenames
lambda = lambda_list;
filenames = {new_filename};

% Subtract out early level response functions
% EEGnew = EEG; % create copy to be modified during run (but only original will be saved)
for iLevel=1:nLevels
    % Get regressor matrices
    if ~exist('Y','var')
        fprintf('%s Getting matrices...',datestr(now,16))
        tic;    
        [X,Y] = GetGlmMatrices(EEG,regressor_events{iLevel},influence{iLevel},artifact_events,artifact_influence{iLevel},demeanX,normalizeX);
        t=toc;
        fprintf('took %.1f seconds.\n',t);
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
    Nh = length(tResponse{iLevel});
    Nr = length(regressor_events{iLevel});
    p = Nh*Nr;
    RF = nan(D,Nh,Nr,numel(lambda));   
    for i=1:numel(lambda)                
        fprintf('%d/%d...\n',i,numel(lambda));        
        beta = (XTX+lambda(i)*eye(p))^(-1)*XTY; % size pxD
        RF(:,:,:,i) = reshape(beta',[D,Nh,Nr]);
    end
    t=toc;
    fprintf('took %.1f seconds.\n',t);
    
    % Place Response functions
    responseFns{iLevel} = RF;       
    
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
    
%     if iLevel<nLevels
%         fprintf('%s Subtracting out results...',datestr(now,16))
%         tic;   
%         Y = Y-RF*X; %TODO: FIX THIS!
%         t=toc;
%         fprintf('took %.1f seconds.\n',t);
%     end
end

% Clean up by deleting lower-level files
fprintf('%s Cleaning up...',datestr(now,16));
tic;
for iLevel = 1:nLevels-1
    if ~strcmp(filenames{iLevel},filenames{nLevels})
        try
            fprintf('Deleting file %s...\n',filenames{iLevel});
            delete(filenames{iLevel});
        catch
            warning('Couldn''t delete file %s!',filenames{iLevel});
        end
    end
end
t=toc;
fprintf('took %.1f seconds.\n',t);

% Update status
fprintf('%s - Analysis complete!\n',datestr(now,16));