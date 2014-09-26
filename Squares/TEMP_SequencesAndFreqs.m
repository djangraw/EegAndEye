% Created 1/21/13 by DJ for one-time use.

subject = 19;
load(sprintf('sq-%d-all.mat',subject)); % y


% nSessions = numel(y);
% [endx endy class] = deal(cell(1,nSessions));
% for i=1:nSessions
%     endx{i} = y(i).saccade.end_position(:,1);
%     endy{i} = y(i).saccade.end_position(:,2);
%     class{i} = y(i).saccade.class;
% %     scatter(endx,endy);
% end
% endx = cat(1,endx{:});
% endy = cat(1,endy{:});
% class = cat(1,class{:});
% 
% Constants = GetSquaresConstants;
% cla; hold on
% % scatter(endx(class==Constants.DISTRACTOR), endy(class==Constants.DISTRACTOR),'g.');
% % isTarget = ismember(class,[Constants.INTEGRATION Constants.COMPLETION]);
% % scatter(endx(isTarget), endy(isTarget),'b.');
% isInt = class==Constants.INTEGRATION;
% isCom = class==Constants.COMPLETION;
% isExt = class==Constants.EXTRA;
% scatter(endx(isInt), endy(isInt),'b.');
% scatter(endx(isCom), endy(isCom),'c.');
% scatter(endx(isExt), endy(isExt),'m.');
% 
% plot(Constants.SQUARE_X, Constants.SQUARE_Y,'ks','MarkerSize',30,'LineWidth',3);
% axis image

%%
% Get sequence of each trial
sequence = cell(1,nSessions);
for i=1:nSessions
    sequence{i} = y(i).trial.is_target_color;
    sequence{i}(y(i).trial.is_right_cross,:) = fliplr(sequence{i}(y(i).trial.is_right_cross,:));
    
end
sequence = cat(1,sequence{:});

% Get unique sequences
seqs = unique(sequence,'rows');
iThisSeq = cell(1,size(seqs,1));
nThisSeq = zeros(1,size(seqs,1));
for j=1:size(seqs,1)
    iThisSeq{j} = find(ismember(sequence,seqs(j,:),'rows'));
    nThisSeq(j) = numel(iThisSeq{j});
end

% sort by prevalence
[nThisSeq,order] = sort(nThisSeq,'descend');
iThisSeq = iThisSeq(order);
seqs = seqs(order,:);

% cla;
% plot(nThisSeq);

%% Get Theta activity on electrode Pz

% Set parameters
chan = 'PZ';
freqrange = [4 9]; % in Hz
epochlims = [-2 5]; % in seconds
colorlims = [0 10]; % in uV?
iChan = strmatch(chan,{EEG.chanlocs.labels}); % get channel number
%% Load data
EEG = pop_loadset(sprintf('sq-%d-all-filtered-50Hz-interpduds.set',subject));

% Filter in theta range
EEG = pop_eegfilt( EEG, freqrange(1), freqrange(2), [], 0); % What kind of filter is this?

% Use Hilbert transform to get power
disp('Getting Hilbert Transform...')
hilb_data = hilbert(EEG.data')';
EEG.data = abs(hilb_data);

% Epoch data
EEG = pop_epoch( EEG, {  'TrialStart-D'  'TrialStart-T'  }, epochlims);

%% Plot result as image

% erpimage(EEG.data(iChan,:,:),[],EEG.times);
figure(210);
clf;
pop_erpimage(EEG,1, iChan,[],chan,10,1,{ 'FixOn-D' 'FixOn-T'},[],'latency' ,'yerplabel','\muV','erp','on','limits',[NaN NaN colorlims NaN NaN NaN NaN],'cbar','on','caxis',colorlims,'topo', { [iChan] EEG.chanlocs EEG.chaninfo } );

%% Plot trial types separately
% get sequence types
iSeq = {find(sum(seqs,2)==0), find(sum(seqs,2)==1), find(sum(seqs,2)==2), find(sum(seqs,2)==3)};
erp = nan(size(EEG.data,1),size(EEG.data,2),numel(iSeq));
for j=1:numel(iSeq)
    % Make temp dataset
    EEG_tmp = pop_selectevent( EEG, 'epoch',cat(1,iThisSeq{iSeq{j}}) ,'deleteevents','off','deleteepochs','on','invertepochs','off');
    % Plot
    figure(210+j); clf;
    pop_erpimage(EEG_tmp,1, iChan,[],chan,10,1,{ 'FixOn-D' 'FixOn-T'},[],'latency' ,'yerplabel','\muV','erp','on','limits',[NaN NaN colorlims NaN NaN NaN NaN],'cbar','on','caxis',colorlims,'topo', { [iChan] EEG.chanlocs EEG.chaninfo } );
    erp(:,:,j) = mean(EEG_tmp.data,3);
    MakeFigureTitle(sprintf('Sequence %d: [%s]',iSeq{j}(1),num2str(seqs(iSeq{j}(1),:))));
end
%%
figure(subject); clf;
plot(EEG.times,squeeze(erp(iChan,:,:))');
hold on
plot([0 0],get(gca,'ylim'),'k')
legend('0 targets','1 target', '2 targets', '3 targets');
xlabel('time from trial start (ms)')
ylabel(sprintf('Power in %d-%d Hz range on electrode %s',freqrange(1),freqrange(2),chan));
title(EEG.setname);
