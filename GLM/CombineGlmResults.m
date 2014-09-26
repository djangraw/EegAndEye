function [combo,responseFns] = CombineGlmResults(results,multcompare)

% Combines the saved results of a GLM analysis from multiple subjects,
% using only the electrodes common to all the results.
%
% combo = CombineGlmResults(results,multcompare)
%
% INPUTS:
% - results is an n-element vector of structs with fields 'responseFns', 
% 'regressor_events', 'tResponse', 'EEG'.
% - multcompare is a string indicating the type of multiple comparisons
% correction you would like to use. ('none','fdr','bonferroni')
%
% OUTPUTS:
% - combo is a struct with the same fields as the input structs and another
% field 'nSubjects' which equals n (the number of inputs), and a field
% 'Pval' which is the p value of each time point relative in a one-sided
% t-test against a mean of zero (with specified multiple comparisons
% correction).
% - responseFns is a 4-D matrix in which each 4th-D page contains the
% response functions from one of the subjects, with only the channels in
% combo.EEG.chanlocs included.
%
% Created 1/4/12 by DJ.
% Updated 7/10/12 by DJ - added p values, multiple comparison corrections,
% and a responseFns output
% Updated 1/22/13 by DJ - use UpdateGlmResultsFomat
% Updated 4/30/13 by DJ - responseFns in cells

if nargin<2
    multcompare = 'fdr';
end

% Get new results format
results = UpdateGlmResultsFormat(results);

% Get list of all channels
channels = {results(1).EEG.chanlocs.labels};
for i=2:numel(results)
    channels = intersect(channels, {results(i).EEG.chanlocs.labels});
end
% Get indices of ok channels
RFcell = cell(1,numel(results));
for i=1:numel(results);
    RFcell{i} = results(i).responseFns{results(i).iLevel}(ismember({results(i).EEG.chanlocs.labels},channels),:,:);
end

% Get p values
disp('Calculating Group Statistics...')
h = waitbar(0,'Calcluating p values...');
responseFns = cat(4,RFcell{:});
Pval_start = nan(size(responseFns,1),size(responseFns,2),size(responseFns,3));
for i=1:size(results(1).responseFns,1)
    waitbar(i/size(results(1).responseFns,1),h);
    for j=1:size(results(1).responseFns,3)
        [~,Pval_start(i,:,j)] = ttest(squeeze(responseFns(i,:,j,:))',0,0.05,'left');
    end
end
close(h);

% Correct for multiple comparisons
switch multcompare
    case 'none'
        disp('No multiple comparisons correction!')
        Pval = Pval_start;
    case 'fdr' % false discovery rate
        disp('Multiple comparisons correction using False Discovery Rate...')
        isHigh = Pval_start>0.5;
        Pval_start(isHigh) = 1-Pval_start(isHigh);        
        Pval = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));
        Pval(isHigh) = 1-Pval(isHigh);
    case 'bonferroni' % bonferroni correction
        disp('Multiple comparisons correction using Bonferroni...')
        isHigh = Pval_start>0.5;
        Pval_start(isHigh) = 1-Pval_start(isHigh);
        Pval = bonf_holm(Pval_start);
        Pval(Pval>0.5) = 0.5;
        Pval(isHigh) = 1-Pval(isHigh);
    otherwise
        error('multiple comparisons method not recognized!');
end

disp('Assembling output...')
% Combine results
combo.regressor_events = results(1).regressor_events;
combo.iLevel = results(1).iLevel;
combo.artifact_events = {'Unknown'};
combo.nuisance_events = {'Unknown'};
% Get average response functions
combo.responseFns = mean(responseFns,4);
% Add common fields to struct
combo.tResponse = results(1).tResponse;
combo.EEG = results(1).EEG;
combo.EEG.data = [];
combo.EEG.chanlocs = results(1).EEG.chanlocs(ismember({results(1).EEG.chanlocs.labels},channels));
% Add combo fields to struct
% combo.nSubjects = numel(results);
combo.Pval = Pval;
disp('Done!')