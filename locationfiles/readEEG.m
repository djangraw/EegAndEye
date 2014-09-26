function [eegdata,events]=readEEG(filename,channel,duration,offset);

%offset should start from zero
%filename='noise.000';
%channel=[1:87];
channel=channel+1;
%duration=36*1024;
%offset=0;

fid=fopen([filename '.dat'],'r','b');
eegdata=[]; events=[];
D=97;  % 96 EEG channels + 1 event

for offseti=1:length(offset),
    %fseek(fid,4*D*offset(offseti)+4*4,'bof');  %before  preprocessing
    fseek(fid,4*D*offset(offseti),'bof');       %after preprocessing

    blklen=2^14;  %=16384
    slices=floor(duration/blklen);
    for i=1:slices,
        slice = fread(fid,[D blklen],'single');
        eegdata=[eegdata slice(channel,:)];
        events=[events slice(1,:)];
        fprintf(['\rReading ' filename '.bin - offset: %d/%d - slice: %d/%d'],offseti,length(offset),i.*blklen,duration);
    end

    blklenlast=duration-(slices.*blklen);
    slice =fread(fid,[D blklenlast],'single');
    eegdata=[eegdata slice(channel,:)];
    events=[events slice(1,:)];

    fprintf(['\rReading ' filename '.bin - offset: %d/%d - slice: %d/%d'],offseti,length(offset),slices.*blklen+blklenlast,duration);

end

fclose(fid);
