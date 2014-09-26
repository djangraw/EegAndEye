function AUC = ClassifierWrapper_forSFFS(alldata,foo,ALLEEG,y,offset)

% Created 8/23/13 by DJ.

useEEG = 0;
bigFeature = alldata';
cvmode = '10fold';

% load(sprintf('ALLEEG%d_NoeogEpochedIcaCropped.mat',subject)); % ALLEEG
% y = loadBehaviorData(subject,sessions,'3DS');

R = ClassifyWithEegAndDwellTime(ALLEEG,y,cvmode,offset,bigFeature,useEEG);

AUC = R.Az;