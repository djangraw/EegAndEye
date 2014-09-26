%% Set up
experiment = 'sq';
switch experiment
    case 'sq'
        basedir = '/Users/dave/Documents/Data/Squares';
        subjects = [9:11, 13:15, 17:19 20:27];
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

%% Load data

%% Recalculate
old_suffix = 'Type-v3pt6-RampUp-ridge014';
new_suffix = 'Type-v3pt6-RampUp-ridgeiter';
for i=1:numel(subjects)
    fprintf('---SUBJECT %d---\n',subjects(i));
    old_filename = sprintf('%s-%d-GLMresults-%s',experiment,subjects(i),old_suffix);
    fprintf('Loading %s...\n',old_filename);
    R = load(old_filename);
    new_filename = sprintf('%s-%d-%s',experiment,subjects(i),new_suffix);
    GetLambdaIteratively(R,2:19,new_filename)
end

