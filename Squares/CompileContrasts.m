function [group_Z, group_RF, legendstr, titlestr] = CompileContrasts(experiments,suffixes,rules,multcorrect)

% Compile various results into a single matrix and apply multiple
% comparsions correction.
%
% [group_Z, group_RF, legendstr] = CompileContrasts(experiments,suffixes,rules,multcorrect)
% 
% INPUTS:
% - experiments is a cell array of strings ('sf','sq','or 'sf3').
% - suffixes is a cell array of strings (e.g., '-Type-v3pt6-RampUp-10fold')
% - rules is a cell array of strings (e.g. 'T0vD0')
%
% OUTPUTS:
% - group_Z is a matrix of Z scores
% - group_RF is a matrix of contrast RFs
% - legendstr is a cell array of strings
%
% Created 9/16/14 by DJ based on GetContrastResults_fast.


if ~exist('multcorrect','var')
    multcorrect = 'none';
end

if ~iscell(experiments)
    experiments = {experiments};
    include_exp = 0;
else
    include_exp = 1;
end
if ~iscell(suffixes)
    suffixes = {suffixes};
    include_suf = 0;
else
    include_suf = 1;
end
if ~iscell(rules)    
    rules = {rules};
    include_rul = 0;
else
    include_rul = 1;
end
%% Reconstruct multiple contrasts

[group_RF,group_Z,group_P,contrastFns] = deal([]);
legendstr = {};

for i=1:numel(experiments)
    experiment = experiments{i};
    for j=1:numel(suffixes)
        suffix = suffixes{j};
        for k=1:numel(rules)
            rule = rules{k};

            if exist(sprintf('%s-%sResults-%scontrast.mat',experiment,suffix,rule),'file')
                R = load(sprintf('%s-%sResults-%scontrast',experiment,suffix,rule));
            else
                warning('DJ:CompileContrasts:NotFound','Results file not found... trying alternate filename.');
                R = load(sprintf('%s-%s-%scontrast',experiment,suffix,rule));
            end
            group_RF = cat(3,group_RF,R.group_RF);
            group_Z = cat(3,group_Z,R.group_Z);
            group_P = cat(3,group_P,R.group_P);
%             contrastFns = cat(3,contrastFns,R.contrastFns);
        %     legendstr = cat(2,legendstr,foo.legendstr);
            newstr = '';
            if include_exp
                newstr = [newstr experiment ', '];
            end
            if include_suf
                newstr = [newstr suffix ', '];
            end
            if include_rul
                newstr = [newstr rule ', '];
            end
            if ~isempty(newstr)
                newstr = newstr(1:end-2);
            end
            legendstr = cat(2,legendstr,newstr);
%             legendstr = cat(2,legendstr,sprintf('%s, %s, %s',experiment,suffix, rule));
        end
    end
end

%% Correct group-level stats
% multcorrect = 'fdr';
clear Pval_end
for i=1:size(group_P,3);
    Pval_start = group_P(:,:,i); 
    if all(isnan(Pval_start(:)))
        Pval_end(:,:,i) = Pval_start;
    else
        switch multcorrect
            case 'fdr' % false discovery rate
                isHigh = Pval_start>0.5;
                Pval_start(isHigh) = 1-Pval_start(isHigh);        
                Pval = reshape(mafdr(Pval_start(:),'bhfdr',true),size(Pval_start));
                Pval(isHigh) = 1-Pval(isHigh);
            case 'none'
                Pval = Pval_start;
        end
        Pval_end(:,:,i) = Pval;
    end
end

group_Z = norminv(Pval_end);

%% Get title string listing consistent elements

titlestr = '';
if ~include_exp
    titlestr = [titlestr experiment ', '];
end
if ~include_suf
    titlestr = [titlestr suffix ', '];
end
if ~include_rul
    titlestr = [titlestr rule ', '];
end
if ~isempty(titlestr)
    titlestr = titlestr(1:end-2);
end