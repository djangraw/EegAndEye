function data_extract(filein,kk,capname)

cd eyetrack
load event_fromeyetrack
    %distractortrialused_EEG       1x525                  16212  cell array
    %fixation_pos                  1x525                  75900  cell array
    %fixation_time                 1x525                  76012  cell array
    %searchstart_end               2x525                   8400  double array
    %trialinfo                     4x525                  16800  double array
    %trialused_EEG               525x1                     4200  double array
cd ..

cd EEG
load event_fromEEG
    %realevents            1x2183219           17465752  double array
    %searchstarttime       1x525                   4200  double array
cd ..

if kk==1;
    channel=2:61;
elseif kk<=3 %an 2cd, and nelson
    channel=1:60;
else          %rewired cap  
    channel=[1 4:60 88 89];
end

duration=1000;

fixation_on_tar{3}=[]; fixation_on_dis{3}=[];

saccade_to_tar{3}=[]; saccade_to_dis{3}=[];

RT_tar{3}=[]; RT_dis{3}=[];

sac_distance_tar_4fix{3}=[];
sac_duration_tar_4fix{3}=[];
sac_distance_dis_4fix{3}=[];
sac_duration_dis_4fix{3}=[];
sac_direction_tar_4fix{3}=[];
sac_direction_dis_4fix{3}=[];

sac_distance_tar_4sac{3}=[];
sac_duration_tar_4sac{3}=[];
sac_distance_dis_4sac{3}=[];
sac_duration_dis_4sac{3}=[];
sac_direction_tar_4sac{3}=[];
sac_direction_dis_4sac{3}=[];


for i=1:size(trialinfo,2)
    if trialused_EEG(i) >0  % is a target trial, and should be used

        dis_sac=sqrt((fixation_pos{i}(2,trialused_EEG(i))-fixation_pos{i}(2,trialused_EEG(i)-1))^2 ...  % distance
                +(fixation_pos{i}(1,trialused_EEG(i))-fixation_pos{i}(1,trialused_EEG(i)-1))^2);
        dur_sac=fixation_time{i}(1, trialused_EEG(i))-fixation_time{i}(2, trialused_EEG(i)-1);          %saccade duration
        
        change_x=fixation_pos{i}(1,trialused_EEG(i))-fixation_pos{i}(1,trialused_EEG(i)-1);
        change_y=fixation_pos{i}(2,trialused_EEG(i))-fixation_pos{i}(2,trialused_EEG(i)-1);
        if change_x>=0 & change_y>=0
                dir_sac=(2*pi-atan(change_y/change_x))*180/pi;
        elseif change_x>=0 & change_y<0
                dir_sac=-atan(change_y/change_x)*180/pi;
        else
                dir_sac=(pi-atan(change_y/change_x))*180/pi;
        end
        
             
        dur=fixation_time{i}(2, end)-fixation_time{i}(1, trialused_EEG(i));       %fixation duration (small saccade negelectted)
        if dur>150
            RT_tar{trialinfo(2,i)}=[RT_tar{trialinfo(2,i)} dur];

            sac_distance_tar_4fix{trialinfo(2,i)}=[sac_distance_tar_4fix{trialinfo(2,i)} dis_sac];
            sac_duration_tar_4fix{trialinfo(2,i)}=[sac_duration_tar_4fix{trialinfo(2,i)} dur_sac];
            sac_direction_tar_4fix{trialinfo(2,i)}=[sac_direction_tar_4fix{trialinfo(2,i)} dir_sac]; 
            
            fixationtime=searchstarttime(i)+fixation_time{i}(1, trialused_EEG(i))-searchstart_end(1,i);
            fixation_on_tar{trialinfo(2,i)}=[fixation_on_tar{trialinfo(2,i)} fixationtime];
        end
        
        dur_pre=fixation_time{i}(2, trialused_EEG(i)-1)-fixation_time{i}(1, trialused_EEG(i)-1);  %previous fixation duration
        if dur_pre>100
            sac_distance_tar_4sac{trialinfo(2,i)}=[sac_distance_tar_4sac{trialinfo(2,i)} dis_sac];
            sac_duration_tar_4sac{trialinfo(2,i)}=[sac_duration_tar_4sac{trialinfo(2,i)} dur_sac];
            sac_direction_tar_4sac{trialinfo(2,i)}=[sac_direction_tar_4sac{trialinfo(2,i)} dir_sac]; 
            
            fixationtime=searchstarttime(i)+fixation_time{i}(2, trialused_EEG(i)-1)-searchstart_end(1,i);  %locked to saccade onset
            saccade_to_tar{trialinfo(2,i)}=[saccade_to_tar{trialinfo(2,i)} fixationtime];
        end
        
   
    elseif length(distractortrialused_EEG{i})>0  % is a distractor trial, and should be used
        for j=1:length(distractortrialused_EEG{i})
                
            dis_sac=sqrt((fixation_pos{i}(2,distractortrialused_EEG{i}(j))-fixation_pos{i}(2,distractortrialused_EEG{i}(j)-1))^2 ...
                        +(fixation_pos{i}(1,distractortrialused_EEG{i}(j))-fixation_pos{i}(1,distractortrialused_EEG{i}(j)-1))^2);
            dur_sac=fixation_time{i}(1, distractortrialused_EEG{i}(j))-fixation_time{i}(2, distractortrialused_EEG{i}(j)-1);
                
            if dur_sac<0
                todisp=sprintf('%d %d\n',i, j);
                disp(todisp);
            end       
            
            change_x=fixation_pos{i}(1,distractortrialused_EEG{i}(j))-fixation_pos{i}(1,distractortrialused_EEG{i}(j)-1);
            change_y=fixation_pos{i}(2,distractortrialused_EEG{i}(j))-fixation_pos{i}(2,distractortrialused_EEG{i}(j)-1);
            if change_x>=0 & change_y>=0
                dir_sac=(2*pi-atan(change_y/change_x))*180/pi;
            elseif change_x>=0 & change_y<0
                dir_sac=-atan(change_y/change_x)*180/pi;
            else
                dir_sac=(pi-atan(change_y/change_x))*180/pi;
            end
            
            
            dur=fixation_time{i}(2, distractortrialused_EEG{i}(j))-fixation_time{i}(1, distractortrialused_EEG{i}(j));
            if dur>150
                RT_dis{trialinfo(2,i)-3}=[RT_dis{trialinfo(2,i)-3} dur];
 
                sac_distance_dis_4fix{trialinfo(2,i)-3}=[sac_distance_dis_4fix{trialinfo(2,i)-3} dis_sac];
                sac_duration_dis_4fix{trialinfo(2,i)-3}=[sac_duration_dis_4fix{trialinfo(2,i)-3} dur_sac];
                sac_direction_dis_4fix{trialinfo(2,i)-3}=[sac_direction_dis_4fix{trialinfo(2,i)-3} dir_sac]; 

                fixationtime=searchstarttime(i)+fixation_time{i}(1, distractortrialused_EEG{i}(j))-searchstart_end(1,i);
                fixation_on_dis{trialinfo(2,i)-3}=[fixation_on_dis{trialinfo(2,i)-3} fixationtime];
            end
            
            dur_pre=fixation_time{i}(2, distractortrialused_EEG{i}(j)-1)-fixation_time{i}(1, distractortrialused_EEG{i}(j)-1);  %previous fixation duration
            if dur_pre>100
                sac_distance_dis_4sac{trialinfo(2,i)-3}=[sac_distance_dis_4sac{trialinfo(2,i)-3} dis_sac];
                sac_duration_dis_4sac{trialinfo(2,i)-3}=[sac_duration_dis_4sac{trialinfo(2,i)-3} dur_sac];
                sac_direction_dis_4sac{trialinfo(2,i)-3}=[sac_direction_dis_4sac{trialinfo(2,i)-3} dir_sac]; 
            
                fixationtime=searchstarttime(i)+fixation_time{i}(2, distractortrialused_EEG{i}(j)-1)-searchstart_end(1,i);  %locked to saccade onset
                saccade_to_dis{trialinfo(2,i)-3}=[saccade_to_dis{trialinfo(2,i)-3} fixationtime];
            end

    
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd mix
for i=1:3
    if kk==1
        eeg_fix2tar{i}=readcogdata(filein,channel,duration,fixation_on_tar{i}-300);
    else
        eeg_fix2tar{i}=readEEG(filein,channel,duration,fixation_on_tar{i}-300);
    end   
    eegtemp=reshape(eeg_fix2tar{i}, length(channel), duration, length(fixation_on_tar{i}));
    [eeg]=pop_importdata('data',eegtemp, 'srate',1000,'pnts',800,'xmin',-0.3, 'nbchan',60,'chanlocs',capname);
    setname=sprintf('eeg-fix2tar-%1d.set', i);
    pop_saveset(eeg,setname);
end

save eeg_fix2tar eeg_fix2tar RT_tar sac_distance_tar_4fix sac_duration_tar_4fix sac_direction_tar_4fix
clear eeg_fix2tar eegtemp eeg

for i=1:3
    if kk==1
        eeg_fix2dis{i}=readcogdata(filein,channel,duration,fixation_on_dis{i}-300);
    else
        eeg_fix2dis{i}=readEEG(filein,channel,duration,fixation_on_dis{i}-300);
    end
    eegtemp=reshape(eeg_fix2dis{i}, length(channel), duration, length(fixation_on_dis{i}));
    [eeg]=pop_importdata('data',eegtemp, 'srate',1000,'pnts',800,'xmin',-0.3, 'nbchan',60,'chanlocs',capname);
    setname=sprintf('eeg-fix2dis-%1d.set', i);
    pop_saveset(eeg,setname);
end
save eeg_fix2dis eeg_fix2dis RT_dis sac_distance_dis_4fix sac_duration_dis_4fix sac_direction_dis_4fix
clear eeg_fix2dis eegtemp eeg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:3
    if kk==1
        eeg_sac2tar{i}=readcogdata(filein,channel,duration,saccade_to_tar{i}-300);
    else
        eeg_sac2tar{i}=readEEG(filein,channel,duration,saccade_to_tar{i}-300);
    end   
    eegtemp=reshape(eeg_sac2tar{i}, length(channel), duration, length(saccade_to_tar{i}));
    [eeg]=pop_importdata('data',eegtemp, 'srate',1000,'pnts',500,'xmin',-0.3, 'nbchan',60,'chanlocs',capname);
    setname=sprintf('eeg-sac2tar-%1d.set', i);
    pop_saveset(eeg,setname);
end

save eeg_sac2tar eeg_sac2tar sac_distance_tar_4sac sac_duration_tar_4sac sac_direction_tar_4sac
clear eeg_sac2tar eegtemp eeg

for i=1:3
    if kk==1
        eeg_sac2dis{i}=readcogdata(filein,channel,duration,saccade_to_dis{i}-300);
    else
        eeg_sac2dis{i}=readEEG(filein,channel,duration,saccade_to_dis{i}-300);
    end
    eegtemp=reshape(eeg_sac2dis{i}, length(channel), duration, length(saccade_to_dis{i}));
    [eeg]=pop_importdata('data',eegtemp, 'srate',1000,'pnts',500,'xmin',-0.3, 'nbchan',60,'chanlocs',capname);
    setname=sprintf('eeg-sac2dis-%1d.set', i);
    pop_saveset(eeg,setname);
end
save eeg_sac2dis eeg_sac2dis sac_distance_dis_4sac sac_duration_dis_4sac sac_direction_dis_4sac


cd ..
