%%% -------------------------------------------------------------------%%%
%%% Sparse Logistic Regression
%%%
%%% Function: Sparse logistic regression
%%% This function performs single trial analysis for EEG data using sparse
%%% logistic regression.
%%%
%%% Input: nothing
%%%     
%%% Output: nothing
%%% 
%%% Jianing Shi
%%% 09/13/2010
%%%
%%% -------------------------------------------------------------------%%%

function LR = SparseLogisticRegression(ALLEEG,setlist,exptname,domain,windowlength,windowoffset,labeltype,chanexclude)

%%% Sparse decoding for channel & time

%%% Add path
addpath(genpath(pwd));

% initializations and checks
 
bootstrap = 0; 
DOKFOLD = 1;
                
if bootstrap, showaz=0; else showaz=1; end;

%setlist = eegb.bisetlist(cls,:);
%exptname = classlabel;
%domain = eegb.domainname;
%windowlength = eegb.windowlength;
%windowoffset = eegb.windowoffset;
%labeltype = eegb.labeltype;

chansubset = 1:min([ALLEEG(setlist(1)).nbchan,ALLEEG(setlist(2)).nbchan]);

chansubset(chanexclude) = [];

vinit=zeros(length(chansubset)+1,1);

show=0;

timedomain = 1; frequencydomain = 0;


%%% ------------ Adjust input from units to indices ----------------

switch domain
    
    case 'time'
        if ~exist('windowlength')||isempty(windowlength), windowlength = 50; end    % default
        if ~exist('windowoffset')||isempty(windowoffset), windowoffset = 0; end     % default
        timedomain = 1; frequencydomain = 0;
        
        fs = ALLEEG(setlist(1)).srate/1000; % samples per ms
        shift = ALLEEG(setlist(1)).xmin*1000; % epoch start relative to stim in ms
        totallength = ALLEEG(setlist(1)).pnts; % samples per epoch
        PSDofEEG = [];
        unit = 'ms';
        
        fprintf(['\n----------------------------------------------------' ...
    '\n   Performing Time Domain Analysis' ...
    '\n----------------------------------------------------\n'])
        
    case 'frequency'
        if ~exist('windowlength')||isempty(windowlength), windowlength = 4; end    % default
        if ~exist('windowoffset')||isempty(windowoffset), windowoffset = 0; end     % default
        timedomain = 0; frequencydomain = 1;
        
        fs = 5;  % samples per Hz                   % JEN: EXPERIMENT HERE
        % specify frequency range for PSD calculation
        frequency = 0:1/fs:100; % in Hz
        shift = frequency(1); % start of PSD relative to 0 freq in Hz
        totallength = length(frequency); % samples in PSD
        PSDofEEG = struct;
        unit = 'Hz';
        
        fprintf(['\n----------------------------------------------------' ...
    '\n   Performing Frequency Domain Analysis' ...
    '\n----------------------------------------------------\n'])

end

% calculate number of trials
ntrial1 = ALLEEG(setlist(1)).trials;
ntrial2 = ALLEEG(setlist(2)).trials;
numtrials = ALLEEG(setlist(1)).trials + ALLEEG(setlist(2)).trials;

% shift back from center to start of window for input to LR plugin
windowstart = windowoffset - windowlength/2; % ms or Hz

% adjust window so it corresponds to indices
windowlength_idx = round(windowlength*fs);
windowstart_idx = round((windowstart-shift)*fs);

% ensure we don't exceed limits of epoch or PSD
windowstart_idx = windowstart_idx(find(windowstart_idx<=totallength-windowlength_idx & windowstart_idx>0));

% subsample windows
numpoints = 10;
windowidx = round(linspace(1,windowlength_idx,numpoints));

% channel number
numchans = length(chansubset);

% window number
numwins = length(windowstart_idx);

%%% ----------------------- make the labels --------------------------
clear bigX1 bigX2;

for widx = 1:length(windowstart_idx)
    
    %LR.windowidx(widx) = widx;
    %LR.windowstart_idx(widx) = windowstart_idx(widx);
    %LR.windowoffset_idx(widx) = windowstart_idx(widx)+round(windowlength_idx/2);
    %LR.windowlength_idx = windowlength_idx;
   
    if timedomain % concatenate the datasets      
         
        x1 = ALLEEG(setlist(1)).data(chansubset,...
             windowstart_idx(widx)-1+windowidx,1:ALLEEG(setlist(1)).trials);
        x2 = ALLEEG(setlist(2)).data(chansubset,...
             windowstart_idx(widx)-1+windowidx,1:ALLEEG(setlist(2)).trials);
        
        x1 = double(x1);
        x2 = double(x2);
        
        x1 = squeeze(mean(x1,2)); x1 = x1';  % [ntrials,numchans]
        x2 = squeeze(mean(x2,2)); x2 = x2';  
         
        y1 = ones(1,ntrial1);
        y2 = -ones(1,ntrial2);
        
        bigX1(:,:,widx) = x1;                                   % [ntrials,numchans,numwins]
        bigX2(:,:,widx) = x2;
             
    elseif frequencydomain % concatenate the PSDs of the datasets
        
        x1 = PSDofEEG(setlist(1)).PSD(chansubset,...
             windowstart_idx(widx)-1+windowidx,1:ALLEEG(setlist(1)).trials);
        
        x2 = PSDofEEG(setlist(2)).PSD(chansubset,...
             windowstart_idx(widx)-1+windowidx,1:ALLEEG(setlist(1)).trials);
         
        x1 = double(x1);
        x2 = double(x2);
        
        x1 = reshape(x1,numchans,numpoints*ntrial1);
        x2 = reshape(x2,numchans,numpoints*ntrial2);
        
        y1 = ones(numpoints,ntrial1);
        y2 = -ones(numpoints,ntrial2);
         
    end
    
end

bigX1 = reshape(bigX1,[ntrial1,numchans*numwins]);           %[ntrials,numchans*numwins]
bigX2 = reshape(bigX2,[ntrial2,numchans*numwins]);

bigX1 = bigX1';                                              %[numchans*numwins,numtrials]
bigX2 = bigX2';

%%% bigX1, bigX2
%%% [numchans*numwins,numtrials]

% Rearrange data for logist.m
% [D x (T x trials)]' OR [D x (F x trials)]

%%% -------------------- Leave One Out Analysis -----------------------
x = cat(2,bigX1,bigX2);
y = cat(2,y1,y2);
clear x1 x2 y1 y2;
clear sol proj;

projloo = zeros(numtrials,1); truthloo = zeros(numtrials,1);
ymean = zeros(numtrials,1); xmean = zeros(numtrials,numchans*numwins);

for loo = 1:numtrials

    % Get leave-out-out data
    idxloo = loo;

    xtrain = x;
    xtrain(:,idxloo) = [];
    xtrain = reshape(xtrain,numchans*numwins,[]);
    xtrain = xtrain';

    ytrain = y; 
    ytrain(:,idxloo) = [];
    ytrain = reshape(ytrain,1,[]);
    ytrain = ytrain';

    xtest = x(:,idxloo);
    xtest = reshape(xtest,numchans*numwins,[]);
    xtest = xtest';

    ytest = y(:,idxloo);
    ytest = reshape(ytest,1,[]);
    ytest = ytest';
    
    %%% Initialization for parameters
    [Opts, Para] = SparseSolverInitialization();

    %%% Sparse logistic regression solver
    Out = SparseLogregSolver(xtrain,ytrain,Opts,Para);
    
    for k = 1:Para.maxStep
       
    % Retrieve the weight
    w = Out.Wflow(:,k);
    v = Out.Vflow(:,k);
    
    % Cardinality
    Cardinality(k) = length(find(w));

    % Get the projection
    proj = xtest*w + v;
    projloo((idxloo-1)+1:idxloo,k) = proj; 

    % Get the truth
    truthloo((idxloo-1)+1:idxloo,k) = (ytest==1);

    % Get mean data for forward model
    ymean(idxloo,k) = mean(proj);
    xmean(idxloo,:,k) = mean(xtest,1);
    
    end % k

end % loo

for k = 1:Para.maxStep
   
    % Get projection from loo's estimate
    bploo = bernoull(1,projloo(:,k));

    % Compute classification accuracy
    Azloo(k) = rocarea(bploo,truthloo(:,k));
     
end
 
[maxAzloo,idxopt] = max(Azloo);
LR.Card = Cardinality(idxopt);

% Azloo
LR.Azloo = max(Azloo);

% Weights
LR.wloo = Out.Wflow(:,idxopt);
LR.vloo = Out.Vflow(:,idxopt);

% Forward model
xmeanw = reshape(xmean(:,:,idxopt),numtrials,numchans,numwins);

for widx = 1:length(windowstart_idx)
   LR.a{widx} = ymean(:,idxopt) \ squeeze(xmeanw(:,:,widx)); 
end

% Print out results
fprintf('Classlabel: %s\t  loo Az: %6.2f\t  card: %d\n',exptname,LR.Azloo,LR.Card);

LR.windowoffset = windowoffset;
LR.chansubset = chansubset;
LR.numchans = length(chansubset);
LR.numwins = numwins;