% clear all;
% close all;
filename = 'pdvo_07_M_070111_15_37pm_03t';
ascfile = 'pdvo_07_M_070111_15_37pm_03t.asc';
[sti_time image_label] = find_events_lb(ascfile,'image');
[time_tmp,eye_xy,pd] = find_events_lb(ascfile,'eyesample');
synctime = find_events_lb(ascfile,'synctime');
% Also record button press time
button_press_times_tmp = find_events_lb(ascfile,'button');

% pd cannot be 0, hence find the first element which has a value of 0, and
% get rid of all the elements afterwards. (rid of all extra elements that
% are initially preallocated)

% find(pd==0,1) -> find the first element which is 0.
% e.g. find(pd==0,2), find the first two elements which are 0
% also, find(pd==0,1,'last')

% get rid of all the 0s in pd, pd(end-rid_index:end) meaning get rid of the
% last rid_index elements of pd. 
rid_index = length(pd)-find(pd==0,1);
pd(end-rid_index:end) = [];

% don't forget to set the length of time to the same as pd.
time_tmp(end-rid_index:end) = [];

% Some of the button press can be subject playing around..hence discard
% those that appeared before the experiment really starts.

iter = 0;

for i = 1:length(button_press_times_tmp)
    if button_press_times_tmp(i) < min(time_tmp)
        button_press_times = button_press_times_tmp(i+1:end);
        iter = iter + 1;
    else
        break;
    end
end

if iter == 0
    % no response is discarded
    button_press_times = button_press_times_tmp;
end

% plot real-time pd
fig1 = figure;

%% interpolate the NaN data points

nan = find(isnan(time_tmp)); 
marker = [];

for i = 1:length(nan)-1
    if nan(i+1)==nan(i)+1
        marker = [marker, nan(i)];
        lower_marker = min(marker)-1;
        upper_marker = max(marker)+1;
        
        if isnan(pd(upper_marker))
            if i~=length(nan)-1
                continue;
            else
                % for the last round, otherwise it would directly follow
                % the continue command, and jump out of the loop
                marker = [marker, nan(i+1)];
                upper_marker = max(marker)+1;
                
                % meaning if the last data element of pd is NaN, then we
                % shall set the last session of NaN to the pd value of the
                % nearest non-NaN data points.
                if upper_marker > length(pd)
                    pd(lower_marker:end) = pd(lower_marker);
                    marker = [];
                else
                    t = lower_marker:(upper_marker-lower_marker):upper_marker;
                    p = pd(lower_marker):(pd(upper_marker)-pd(lower_marker)):pd(upper_marker);
                    x = lower_marker:upper_marker;
                    pd(lower_marker:upper_marker) = interp1(t,p,x);
                end
                break;
            end
        end
        
    else
        if lower_marker==0
            % if NaN is from the very beginning, which means lower_marker
            % would be 0.
            marker = [marker,nan(i)];
            upper_marker = max(marker)+1;
            pd(1:upper_marker) = pd(upper_marker);
            marker = [];
        else
            marker = [marker,nan(i)];
            lower_marker = min(marker)-1;
            upper_marker = max(marker)+1;
            
            t = lower_marker:(upper_marker-lower_marker):upper_marker;
            
            % if accidentally pd(upper_marker) = pd(lower_marker), then
            % simply set p as below. For under such circumstances, the
            % difference would be 0, and would end up in an empty matrix
            % (i.e. 2:0:2 -> [], 2:2 cannot work, since it would only
            % produce 1 element, hence simply put the two back together.
            if pd(upper_marker) == pd(lower_marker)
                p = [pd(lower_marker),pd(upper_marker)];
            else
                p = pd(lower_marker):(pd(upper_marker)-pd(lower_marker)):pd(upper_marker);
            end
            x = lower_marker:upper_marker;
            pd(lower_marker:upper_marker) = interp1(t,p,x);
            marker = [];      
        end
    end
end

time = min(time_tmp) + (0:(length(pd)-1));
plot(time-min(time),pd/mean(pd));
hold on;
plot(nan,pd(nan)/mean(pd),'mx');

% % plot events
hold on;