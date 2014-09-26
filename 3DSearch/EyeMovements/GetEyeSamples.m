function GetEyeSamples(subject,sessions,prefix,underscore)

% Saves the raw eye positions and pupil sizes so we don't have to re-run
% them later.
%
% GetEyeSamples(subject,sessions,prefix,underscore)
%
% Run this program before running PlotEyeErps_MultiSession.
%
% INPUTS:
% -subject is a scalar indicating the number of the subject being recorded.
% -sessions is a vector of session numbers.
% -prefix is a string: your filename before the first hyphen. (default:
% '3DS')
% -underscore is a binary value that is 1 if the hyphens in the normal 
% input filename are replaced with underscores. (default: 0)
%
% The path must have access to files named 
% '<prefix>-<subject>-<sessions(i)>.asc' containing the raw eye position 
% and pupil size data.
% This program will save files named 
% '<prefix>-<subject>-<sessions(i)>-eyepos.mat' to the current directory.
%
% Created 6/3/11 by DJ.
% Updated 2/10/12 by DJ - added prefix/underscore inputs for squares data

% Set up
if nargin<3
    prefix = '3DS';
end
if nargin<4
    underscore = 0;
end
    
% Import and save data
for i=1:numel(sessions)    
    fprintf('Importing file %d of %d...\n',i,numel(sessions));
    % Get filenames
    if underscore
        filename_in = sprintf('%s_%d_%d.asc',prefix,subject,sessions(i));
    else
        filename_in = sprintf('%s-%d-%d.asc',prefix,subject,sessions(i));
    end
    filename_out = sprintf('%s-%d-%d-eyepos.mat',prefix,subject,sessions(i));
    % Import data
    [eyepos, pupilsize] = find_events(filename_in,'eyesample');
    % Save results
    save(filename_out,'eyepos','pupilsize');
end

disp('Done!')