experiment = 'sf3';
subjects = 1:12;

% experiment = 'sf';
% subjects = [1:10 12:13];

% experiment = 'sq';
% subjects = [9:11, 13:15 17:27];

% old_suffix = 'GLMResults-Type-v3pt6-RampUp';
old_suffix = 'GLMResults-Type-v3pt6-Peak';
% new_suffix = 'Type-v3pt6-RampUp-Matrices';
% new_suffix = 'Type-v3pt6-Matrices';
new_suffix = 'Type-v3pt6-Peak-Matrices';

if strcmp(experiment,'sf3')
    iEvents = 2:25;
%     iEvents = 2:14;
else
    iEvents = 2:19;
%     iEvents = 2:12;
end
demeanX = 0;
demeanY = 0;
normalizeX = 1;
normalizeY = 1;


switch experiment
    case 'sq'
        basedir = '/Users/dave/Documents/Data/Squares';
        subjects = [9:11, 13:15, 17:27];
        folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
            '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
            '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
            '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
            '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
            '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};
    case 'sf'
        basedir = '/Users/dave/Documents/Data/SquaresFix';
        subjects = [1:10 12:13];
        folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
            '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
            '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
            '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
    case 'sf3'
        basedir = '/Users/dave/Documents/Data/SquaresFix3';
        subjects = 1:12;
        folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
            '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
            '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};
end


%% GET MATRICES AND SAVE THEM

for iSubj=1:numel(subjects)
    fprintf('---%s Subject %d---\n',experiment, subjects(iSubj));
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
    [X,Y,Xmean,Xrss] = GetGlmMatrices(R.EEG,R.regressor_events{1}(iEvents),R.influence{1},R.artifact_events,R.artifact_influence{1},demeanX,normalizeX);
    
    % Normalize Y
    fprintf('Normalizing Y...\n');
%     Ymean = mean(Y,1); 
    Ymean = zeros(1,size(Y,2)); % don't de-mean Y~
    Yrss = nan(1,size(Y,2));
    for i=1:size(Y,2)
        foo = Y(:,i)-Ymean(i);
        Yrss(i) = sqrt(foo'*foo);
        Y(:,i) = foo/Yrss(i);
    end
    
    if ~demeanX
        X = sparse(X); % to save space        
    end
    
    %% Save
    fprintf('%s Saving Matrices...\n',datestr(now,16));
    save(sprintf('%s-%d-%s.mat',experiment,subjects(iSubj),new_suffix),'X','Y','Xmean','Ymean','Xrss','Yrss');
    disp('Done!')
end