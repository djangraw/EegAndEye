% function ConsolidateGlmResults_script

% prefix = 'sq';
% subjects = [9:11 13:15 17:19];
prefix = 'sf';
subjects = 1:10;
bands = {'','Theta_v2-','UpperAlpha_v2-','Lower1Alpha_v2-','Lower2Alpha_v2-'};
contrast_suffix_cell = {'','CvsD1' 'CvsI','IvsD0','CvsD0','CvsInc','IncvsD1','D1vsD0'};
switch prefix
    case 'sq'
        old_suffixes_cell = {{'LREvents-v2pt3','LREvents-SqNum-v2pt3','LREvents-SqNum-NewType-v2pt3'},...
            {'LREvents-v2pt3','LREvents-NewType-v2pt3','LREvents-NewType-SqNum-v2pt3'}};
        new_suffix_cell = {'LREvents-SqNum-Type-v2pt4',...
            'LREvents-Type-SqNum-v2pt4'};
    case 'sf'
        old_suffixes_cell = {{'SqNum-v2pt3','SqNum-Type-v2pt3'},...
            {'Type-v2pt3','Type-SqNum-v2pt3'}};
        new_suffix_cell = {'SqNum-Type-v2pt4',...
            'Type-SqNum-v2pt4'};
end
%% Consolidate
for j=1:numel(bands)
    for k=1:numel(contrast_suffix_cell)
        for i=1:numel(new_suffix_cell)
            % set up suffixes
            old_suffixes = old_suffixes_cell{i};
            new_suffix = new_suffix_cell{i};
            for l=1:numel(old_suffixes)
                old_suffixes{l} = strcat(bands{j},old_suffixes{l},contrast_suffix_cell{k});
            end
            new_suffix = strcat(bands{j},new_suffix,contrast_suffix_cell{k});

            % Check if this consolidation has already been done, or if
            % the completed results aren't there to consolidate
            new_filename = sprintf('%s-%d-GLMresults-%s.mat',prefix,subjects(end),new_suffix);
            old_filename = sprintf('%s-%d-GLMresults-%s.mat',prefix,subjects(end),old_suffixes{end});
            if exist(new_filename,'file') || ~exist(old_filename,'file')
                fprintf('--- SKIPPING %s ---\n',new_suffix);
            else % if not, run it!
                fprintf('--- RUNNING %s ---\n',new_suffix);
                ConsolidateGlmResults(prefix,subjects,old_suffixes,new_suffix);
            end
        end
    end
end
%% Remove old