
% Created 9/15/14 by DJ.
% Updated 1/28/15 by DJ - Made event selection more automatic, added comments.
% Updtaed 2/25/15 by DJ - added events output and made sqfix 'squares' rule 1:3

% experiment = 'sf3';
experiment = 'sf';
% experiment = 'sq';

% old_suffix = 'GLMResults-Type-v3pt6';
% old_suffix = 'GLMResults-Type-v3pt6-RampUp';
% old_suffix = 'GLMResults-Type-v3pt6-Peak';
old_suffix = 'GLMResults-SqNum-v3pt6';
% old_suffix = 'GLMResults-TargDis-v3pt6';

% new_suffix = 'Type-v3pt6-Matrices'; %_500ms
% new_suffix = 'Type-v3pt6-RampUp-Matrices';
% new_suffix = 'Type-v3pt6-Peak-Matrices';
new_suffix = 'SqNum-v3pt6-Matrices';
% new_suffix = 'TargDis-v3pt6-Matrices'; %_500ms
% new_suffix = 'Square-v3pt6-Matrices';

if strcmp(new_suffix(1:5),'SqNum')
    iEvents = 1:7; % for SqNum
elseif strcmp(new_suffix(1:7), 'TargDis')
    iEvents = 1:4; % for TargDis
elseif strcmp(new_suffix(1:6),'Square')
%     if strcmp(experiment,'sq')
        iEvents = 1:3; % include TrialStart
%     else
%         iEvents = 1:2; % TrialStart coincides with 1st square event
%     end
elseif strcmp(experiment,'sf3') % sf3/type
    if strcmp(new_suffix(end-12:end-9),'Peak') || strcmp(new_suffix(end-14:end-9),'RampUp')
        iEvents = 2:25; % for peak/rampup
    else
        iEvents = 2:14; % for type
    end
else % sq/type or sf/type
    if strcmp(new_suffix(end-12:end-9),'Peak') || strcmp(new_suffix(end-14:end-9),'RampUp')
        iEvents = 2:19; % for peak/rampup
    else
        iEvents = 2:12; % for type
    end
end

% declare params
demeanX = 0;
demeanY = 0;
normalizeX = 0;%1;
normalizeY = 0;%1;
influence = [0 1000];
% influence = [0 500];

% Get constants
[subjects,basedir,folders] = GetSquaresSubjects(experiment);


%% GET MATRICES AND SAVE THEM

for iSubj=1:numel(subjects)
    fprintf('---%s Subject %d/%d---\n',experiment, iSubj, numel(subjects));
    % Check if we've done this subject already.
    cd(basedir);
    cd(folders{iSubj});
    if exist(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),new_suffix),'file')
        fprintf('Already Done... Skipping!\n');
        continue;
    end
    
    fprintf('%s Subject %d: Loading...\n',datestr(now,16),subjects(iSubj))
    R = load(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),old_suffix));
    
    fprintf('%s Getting Matrices...\n',datestr(now,16))
    [X,Y,Xmean,Xrss] = GetGlmMatrices(R.EEG,R.regressor_events{1}(iEvents),influence,R.artifact_events,R.artifact_influence{1},demeanX,normalizeX);
    
    % Normalize Y
    
    if demeanY
        fprintf('Demeaning Y...\n')
        Ymean = mean(Y,1); 
    else
        Ymean = zeros(1,size(Y,2)); % don't de-mean Y~
    end
    if normalizeY
        fprintf('Normalizing Y...\n');        
    end
    Yrss = ones(1,size(Y,2));
    for i=1:size(Y,2)
        foo = Y(:,i)-Ymean(i);
        if normalizeY
            Yrss(i) = sqrt(foo'*foo);
        end
        Y(:,i) = foo/Yrss(i);
    end
    
    
    if ~demeanX
        X = sparse(X); % to save space        
    end
    
    events = R.regressor_events{1}(iEvents);
    %% Save
    fprintf('%s Saving Matrices...\n',datestr(now,16));
    save(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),new_suffix),'X','Y','Xmean','Ymean','Xrss','Yrss','events');
    disp('Done!')
end