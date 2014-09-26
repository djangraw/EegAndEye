function [fwdmodel time Azloo] = MakeLogRegMovie(ALLEEG,setlist,chansubset,trainingwindowlength,trainingwindowinterval,LOO,bootstrap)

% Uses MakeTopoMovie specifically for the results of pop_logisticregression.
%
% [fwdmodel time looAz] = MakeLogRegMovie(ALLEEG,setlist,chansubset,trainingwindowlength,trainingwindowinterval,LOO)
%
% Inputs (this list taken from pop_logisticregression's help file):
% Inputs:
%   ALLEEG               - array of datasets
%   setlist              - list of datasets (of size 2)
%   chansubset           - vector of channel subset for datasets           
%   trainingwindowlength - Length of training window in samples            
%   trainingwindowinterval - shift of training window(s) in samples
%   LOO                  - do Leave-One-Out analysis
% Inputs to pop_logisticregression hard-coded below:
%   regularize           - regularize [1|0] -> [yes|no]                    [1]
%   lambda               - regularization constant for weight decay.       [1e-6]
%                            Makes logistic regression into a support 
%                            vector machine for large lambda
%                            (cf. Clay Spence)
%   lambdasearch         - [1|0] -> search for optimal regularization      [1]
%                            constant lambda
%   eigvalratio          - if the data does not fill D-dimensional         [1e-4]       
%                            space, i.e. rank(x)<D, you should specify 
%                            a minimum eigenvalue ratio relative to the 
%                            largest eigenvalue of the SVD.  All dimensions 
%                            with smaller eigenvalues will be eliminated 
%                            prior to the discrimination. 
%   vinit                - initialization for faster convergence           [zeros(D+1,1)]
%   show                 - display discrimination boundary during training [0]
% Outputs:
%   fwdmodel - forward model of logistic regression, defined as a=y\x.
%   time  - Vector of epoch time indices
%   Azloo - Vector of LOO Az value for each time
%
% See pop_logisticregression for more info on how the logistic regression is done.
% See MakeTopoMovie for more info on how the data is plotted.
%
% Created 8/12/10 by DJ.
% Updated 8/17/10 by DJ - added LOO analysis (from pop_logisticregression).
% Updated 8/18/10 by DJ - comments
% Updated 9/30/10 by DJ - made LOO optional with optional input LOO
% Updated 11/3/10 by DJ - added Azloo output and comments
% Updated 11/23/10 by DJ - fixed LOO=0 bug
% Updated 1/12/11 by DJ - added bootstrap as input

%% Set up
if nargin<6
    LOO = 1; % do LOO analysis by default
end
if nargin <7
    bootstrap=0;
end
EEG = ALLEEG(setlist(1));
trainingwindowoffset = 1 : trainingwindowinterval : EEG.pnts-trainingwindowlength;
nWindows = numel(trainingwindowoffset);
iMidTimes = trainingwindowoffset + floor(trainingwindowlength/2); % middle of time window
time = EEG.times(iMidTimes)*0.001; % crop and convert to seconds

% Set parameters
regularize = 1;
lambda = 1e-6;
lambdasearch = true;
eigvalratio = 1e-4;
vinit = zeros(length(chansubset)+1,1);
show = 0;

%% Perform logistic regression (results will be saved in ICA data)
if bootstrap
    disp('---Skipping logistic regression...');
    fwdmodel = [];
else
    disp('---Performing logistic regression...');
    ALLEEG = pop_logisticregression(ALLEEG,setlist,chansubset,chansubset,trainingwindowlength,trainingwindowoffset,regularize,lambda,lambdasearch,eigvalratio,vinit,show,0);
    fwdmodel = ALLEEG(setlist(1)).icawinv; % forward model, defined in pop_logisticregression as a=y\x.
end

if LOO
    %% Perform leave-one-out analysis
    disp('---Performing LOO analysis...');
    tic
    N=ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials;
    truth=[zeros(trainingwindowlength.*ALLEEG(setlist(1)).trials,1); ones(trainingwindowlength.*ALLEEG(setlist(2)).trials,1)];
        
    if bootstrap
        disp('BOOTSTRAPPING!!!!!')
        truth=truth(randperm(numel(truth))); % Scramble labels (should we be choosig random trials with replacement?)
    end
    
    for i=1:length(trainingwindowoffset),
        x=cat(3,ALLEEG(setlist(1)).data(chansubset,trainingwindowoffset(i):trainingwindowoffset(i)+trainingwindowlength-1,:), ...
            ALLEEG(setlist(2)).data(chansubset,trainingwindowoffset(i):trainingwindowoffset(i)+trainingwindowlength-1,:));
        x=x(:,:)'; % Rearrange data for logist.m [D (T x trials)]'
    %     ploo=zeros(N*trainingwindowlength,1);

        parfor looi=1:length(x)/trainingwindowlength,
            indx=ones(N*trainingwindowlength,1);
            indx((looi-1)*trainingwindowlength+1:looi*trainingwindowlength)=0;
            tmp = x(indx==1,:); % LOO data
            tmpt = [truth(indx==1)];   % Target

            vloo(:,looi)=logist(tmp,tmpt,vinit,show,regularize,lambda,lambdasearch,eigvalratio);
            y(:,looi) = [x((looi-1)*trainingwindowlength+1:looi*trainingwindowlength,:) ones(trainingwindowlength,1)]*vloo(:,looi);

            ymean(looi)=mean(y(:,looi));

            %      ploo(find(1-indx)) = bernoull(1,y(:,looi));

    %         ploo((looi-1)*trainingwindowlength+1:looi*trainingwindowlength) = bernoull(1,y(:,looi));
            ploomean(looi)=bernoull(1,ymean(looi));

        end
        truthmean=([zeros(ALLEEG(setlist(1)).trials,1); ones(ALLEEG(setlist(2)).trials,1)]);
    %     [Azloo,Ryloo,Rxloo] = rocarea(ploo,truth);
        [Azloo(i),Ryloo,Rxloo] = rocarea(ploomean,truthmean);
        fprintf('Window Onset: %d; LOO Az: %6.2f\n',trainingwindowoffset(i),Azloo(i));
        if Azloo(i)==0
            disp('weird.')
        end
    end

    %% Plot LOO Results
    if ~bootstrap
        figure; hold on;
        plot(time,Azloo);
        plot(get(gca,'XLim'),[0.5 0.5],'k--');
        plot(get(gca,'XLim'),[0.75 0.75],'k:');
        ylim([0.3 1]);
        title('Leave-one-out analysis');
        xlabel('time (s)');
        ylabel('LOO Az');
    end
    toc
else
    Azloo = [];
end

%% Make TopoMovie figure
weights = ALLEEG(setlist(1)).icaweights';
h = MakeTopoMovie(weights,time,EEG.chanlocs);
ylabel(h.ERP,'Spatial Weights');

% h = MakeTopoMovie(fwdmodel,time,EEG.chanlocs);
% ylabel(h.ERP,'Fwd Model');

