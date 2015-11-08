function [subjects,basedir,folders,offsets] = GetSquaresSubjects(experiment)
% [subjects,basedir,folders,offsets] = GetSquaresSubjects(experiment)
%
% INPUTS:
% - experiment is a string ('sq','sf',or 'sf3') indicating the type of
% experiment.
% 
% OUTPUTS:
% - subjects is an N-element vector of subject numbers.
% - basedir is a string indicating the base directory in which the
% subjects' data folders will be found.
% - folders is an N-element cell array of strings indicating each subject's
% data folder name. So the data for subject <subjects(i)> will be in the
% directory <basedir>/folders{i}.
% - offsets is an N-element vector of doubles indicating the offset (in ms)
% of the EEG from the eyelink events (as determined previously by 
% GetTimingOffsets.m).
%
% Created 9/16/14 by DJ.
% Updated 1/15/15 by DJ - fit directories on new NIH machine.
% Updated 1/27/15 by DJ - added offsets output.
% Updated 2/19/15 by DJ - allow long-form (e.g. 'Squares') experiment names

switch experiment
    case {'sq','Squares'}
%         basedir = '/Users/dave/Documents/Data/Squares';
        basedir = '/Users/jangrawdc/Documents/LiincData/Squares';
        subjects = [9:11, 13:15, 17:27];
        folders = {'2012-03-29-Pilot-S9', '2012-04-09-Pilot-S10', '2012-04-18-Pilot-S11', ...
            '2012-05-04-Pilot-S13', '2012-05-29-Pilot-S14', '2012-06-12-Pilot-S15', ...
            '2012-06-15-Pilot-S17', '2012-06-21-Pilot-S18', '2012-06-22-Pilot-S19',...
            '2012-06-25-Pilot-S20', '2013-04-15-Pilot-S21', '2013-04-16-Pilot-S22', ...
            '2014-02-13-Pilot-S23', '2014-02-17-Pilot-S24', '2014-02-20-Pilot-S25', ...
            '2014-03-04-Pilot-S26', '2014-03-10-Pilot-S27'};
%         offsets = [42 102 -12 -18 -12 -12 6 -12 18 0 0 0 0 0 0 0 0];
        offsets = [30 100 -10 -30 -10 -30 0 0 0 -10 0 40 -20 -20 -10  10  20];

    case {'sf','SquaresFix'}
%         basedir = '/Users/dave/Documents/Data/SquaresFix';
        basedir = '/Users/jangrawdc/Documents/LiincData/SquaresFix';
        subjects = [1:10 12:13];
        folders = {'2013-03-05-Pilot-S1','2013-03-11-Pilot-S2','2013-03-18-Pilot-S3',...
            '2013-03-19-Pilot-S4','2013-04-01-Pilot-S5','2013-04-02-Pilot-S6',...
            '2013-04-09-Pilot-S7','2013-04-11-Pilot-S8','2013-04-11-Pilot-S9',...
            '2013-04-17-Pilot-S10','2014-02-17-Pilot-S12','2014-03-06-Pilot-S13'};
%         offsets = zeros(1,numel(subjects));
        offsets = [40 30 70 70 20 0 20 -20 -20 -10 -10 0];

    case {'sf3','SquaresFix3'}
%         basedir = '/Users/dave/Documents/Data/SquaresFix3';
        basedir = '/Users/jangrawdc/Documents/LiincData/SquaresFix3';
        subjects = 1:12;
        folders = {'2013-05-15-Pilot-S1'    '2013-12-02-Pilot-S2'    '2013-12-04-Pilot-S3'    '2013-12-04-Pilot-S4',...
            '2013-12-05-Pilot-S5'    '2013-12-16-Pilot-S6'    '2013-12-19-Pilot-S7'    '2014-02-04-Pilot-S8',...
            '2014-02-06-Pilot-S9'    '2014-02-06-Pilot-S10'   '2014-02-12-Pilot-S11'   '2014-02-18-Pilot-S12'};
%         offsets = zeros(1,numel(subjects));
        offsets = [-20 50 90 80 80 -10 0 -30 20 0 30 20]; % added 7/21/14

end