% Significance test for logisitic regression results.
%
%
%
% Created 1/12/11 by DJ.

%% Plot
doLOO = true; % perform leave-one-out analysis
bootstrap = true; % perform bootstrapping
samplesPerWindow = 13; % at sampling rate of 256: 50.8ms
samplesPerShift = 3; % at sampling rate of 256: 11.7ms

nIter = 100;

for iSubject = 1:3    
    for iter = 1:nIter
        fprintf('===Subject %d, Iteration %d===\n', iSubject, iter);
        
        whichDatasets = (iSubject-1)*2 + [1 2]; % KLUGE       
        [~, time, Azloo] = MakeLogRegMovie(ALLEEG,whichDatasets,1:ALLEEG(whichDatasets(1)).nbchan,samplesPerWindow,samplesPerShift,doLOO,bootstrap);

        % Save the results in the 'LOO.mat' file
        SaveLooResults(ALLEEG(whichDatasets(1)),ALLEEG(whichDatasets(2)),time,Azloo,samplesPerWindow,samplesPerShift,datestr(now));
    end
end