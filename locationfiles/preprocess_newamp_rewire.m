%function preprocess(dirname)
  
% preprocess(dirname) will generate a new data file with: a) clean event
% channel discretized between 1..256. b) subtracts eye blink, horizontal,
% and vertical eye motion activity. c) removes low frequency drifts, 60Hz
% noise, and 120Hz harmonic noise. e) down-samples the data if required. d)
% one can specify the channels that are to be filtered. The rest is left
% alone except for potentially down-sampling. e) the output data file can
% include optionally additional 3 channels with eye activity (after event
% channel). f) The filters can be applied as a FIR with linear phase
% (constant delay). In that case the events channel and unaffected channels
% are delayed by the appropriate amount. This guarantees zero phase
% delay. "dirname" specifies the directory that contain files with eye blink
% calibration data, event calibration data, and the raw data to be
% filtered. Modify the code to set the corresponding file names. Similarly,
% the options outlines above can be set by editing the corresponding
% parameters in the code.
		
% input and output filenames
fileeye   = 'eyecalibrate';
filein    = 'experiment';
fileout   = 'experiment_filtered';
append = 0;  % =1 if you want to append to existing output file		    

% parameters 
fsref = 1000;   % sampling frequency for preocessing and output
show = 1;       % set either 0 or 1 to see stuff
filter60Hz = 1; % select 60Hz notch filtering if needed
filter120Hz = 1; % select 120Hz harmonic notch filtering if needed
filterDC = 1;   % high pass filter to remove DC
blklen = 2^14;  % processing block length - pick as large as memory allows
FIR  = 1;       % =1 use linear phase FIR and delay events by the group delay.
                % =0 use IIR and do not dealy
saveEyeCoord=0; % =1 will save extra eye coordinate as the last 3
                % coordiantes. Note that event channels no longer in last
                % channel, and that file larger than original. =0 simply
		% subtract and leave data otherwise in original format.

% pick the channels to be filtered
channels   = [1 4:89]+1;
%channels   = [1:87]+1;

D=97;
fs=1000;
fileinfo=dir([filein '.dat']);
N=(fileinfo.bytes-16)/97/4;

% sampling frequency ratio relative to fsref
%[D,N,fs,gain] = readheader(filein);
%
% channels not in "channels" variable will be left alone
channelsNo = 1:D; channelsNo(channels) = 0; 
channelsNo = channelsNo(find(channelsNo));

% ratio of given and desired fs
fsr = fs/fsref;

% decimation filter (same as in decimate())
if fsr==1
  adec=1;bdec=1;
else
  [bdec,adec]=cheby1(4,0.05,.8/fsr);
end

% Notch (2nd order Butterworth, bandstop, F3db1=58, F3db2=62)
if filter60Hz 
  [notchnum,notchdenom]=butter(2,[58 62]/fs*2,'stop');
else
  notchnum=1;notchdenom=1;
end

% Notch (2nd order Butterworth, bandstop, F3db1=118, F3db2=122)
if filter120Hz 
  [notchnum120,notchdenom120]=butter(2,[118 122]/fs*2,'stop');
else
  notchnum120=1;notchdenom120=1;
end

% High-pass filter (2nd order Butterworth, cutoof f = 0.5 Hz, fs = 250 Hz).
if filterDC
  [hpnum,hpdenom]=butter(2,0.5/fs*2,'high');
else
  hpnum=1;hpdenom=1;
end

% combine all three filters
a = poly([roots(adec);roots(hpdenom);roots(notchdenom);roots(notchdenom120)]);
b = conv(notchnum120,conv(notchnum,conv(hpnum,bdec)));
delay = 1; % this is how much the event channel is to be delayed

% make linear phase FIR filter from this IIR
if FIR
  u = zeros(2*1024,1); u(round(length(u)/2))=1;
  b=filtfilt(b,a,u); a=1;
  delay = [zeros(1,floor(length(b)/2)-1) 1];
end

if show % the filters we will apply
  freqz(b,a,1024*4,fs); 
  subplot(2,1,2); grpdelay(b,a,1024,fs); 
  if ~FIR, axis([0 fs/2 0 fs/50]); end
  drawnow
end

%keyboard
%tmp=randn(2000,1);tmp(451:550)=randn(100,1)+5; plot([filter(delay,1,tmp)  filter(b,a,tmp)])

% calibrate event markers and eye blink and motions
%[m,s]= calibrate_aim4an(fileevent,100,100); 

fileinfo=dir([fileeye '.dat']);
Neye=(fileinfo.bytes-16)/97/4;
[eyedata, eyeevents]=readEEG_b4preprocessing(fileeye,channels-1,Neye,0) ;
usedcha=[1:size(eyedata,1)];
V(channels,:) = eyecalibrate_newamp(eyedata,eyeevents,usedcha); drawnow

% these channels are not to be included in the projections.
V(channelsNo,:) = 0; 

% this will project out eye artifacts 
U = (eye(size(V,1))-V*pinv(V)); 

% ... and put them in separate coordinates (will dump later if not needed)
U = [U; pinv(V)];

% open i/o files and set counter
fidin  = fopen([filein '.dat'],'r','b');
fidout = fopen([fileout '.dat'],'w','b');
count = 0;

% count-up the block lenghts we have read
Ntotal = count+round(N/fsr);

% Init memory for power tracking in eye blink removal 
W = []; P = []; Z = []; Zdec = []; Zdel = [];

fseek(fidin, 4*4,'bof');

for i=1:ceil(N/blklen)
  
  % read in data 
  X = fread(fidin,[D blklen],'single');

  % filter some delay other and downsample
  [X(channels,:),  Zdec] = filter(    b, a, X(channels,:),   Zdec,2); 
  [X(channelsNo,:),Zdel] = filter(delay, 1, X(channelsNo,:), Zdel, 2);
  X = X(:,fsr:fsr:end);  

  % remove eye blink and motion and add eye coordinates
  X = U*X; 

  % stuff to keep the user informed
  count = count + size(X,2);
  show = (show>0)*size(X,2); 
  disp([num2str(i) ' ' num2str(Ntotal-count)]);
  
  if ~saveEyeCoord, X = X(1:D,:); end;
  % save all data (at new sampling rate)
  fwrite(fidout, X, 'single');
  %if max(max(abs(X)))>=2^15, error('Out of bounds 16bit'); end
  
end

fclose(fidin);
fclose(fidout);













