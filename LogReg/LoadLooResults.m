function LooOut = LoadLooResults(varargin)

% Gets the results of a leave-one-out cross-validation from the file
% 'LOO.mat'.
%
% fitsAll = LoadLooResults(varargin);
%
% POSSIBLE INPUTS:
% -setname1 and setname2 are the names of the datasets e.g.
% '3DS-2-all-targapp'
% -time is a vector of times at the center of the training windows
% -Az is a vector of the leave-one-out classification Az values for each
% window.
% -trainingwindowlength is the length of each training window (in samples).
% -trainingwindowinterval is the distance between each training window
% center. [unique(diff(time))]
% -datetime is the date and time at which the results were saved. [datestr(now)]
% -reference is a cell array containing strings specifying which channels
% were used as a reference. [EEG.ref]  Or the number of references used.
% -filename is the name of the LOO file.
%
% Created 1/12/11 by DJ
% Updated 1/14/11 by DJ - added filename input
% Updated 2/17/11 by DJ - added number-of-references

%% SET UP
% handle inputs
filename = 'LOO.mat'; % default
iField = 0;
for i=1:nargin/2
    if strcmp(varargin{2*i-1},'filename')
        filename = varargin{2*i};
    else
        iField = iField + 1;
        name{iField} = varargin{2*i-1};
        value{iField} = varargin{2*i};
    end
end

% load loo results
fprintf('Getting LOO info from log %s\n',filename);
looLog = which(filename); % finds 'LOO.mat' in the current path
load(looLog); % loads the variable 'LOO'

%% FIND INFO
fitsAll = ones(1,numel(LOO));
for i=1:numel(name)    
    fits = ones(1,numel(LOO));
    switch name{i}
        case {'setname1', 'setname2'}
            % Does the setname contain the search string anywhere in it?
            vals = {LOO(:).(name{i})};
            fits = ~cellfun('isempty',strfind(vals,value{i}));
       
        case {'reference' 'time' 'Az'}
            % Does the array match the search value exactly? ...or does the
            % number of elements in the array match the search value?
            vals = {LOO(:).(name{i})};
            for j=1:numel(vals)
                if isequal(vals{j},value{i})
                    fits(j)=1;
                elseif isequal(length(vals{j}),value{i})
                    fits(j)=1;
                else
                    fits(j)=0;
                end
            end
        case {'trainingwindowlength' 'trainingwindowinterval'}
            % Does the value match the search value exactly?
            vals = [LOO(:).(name{i})];
            fits = (vals==value{i});
        case 'datetime'
            % Is the datetime within 1 day of the given datetime?
            vals = {LOO(:).(name{i})};
            for j=1:numel(vals)
                if abs(datenum(vals{j}) - datenum(value{i})) < 1
                    fits(j)=1;
                else
                    fits(j)=0;
                end
            end
        otherwise
            warning('Field name %s not recognized!',name{i});         
    end
    fitsAll = fitsAll & fits;
end

LooOut = LOO(fitsAll);


