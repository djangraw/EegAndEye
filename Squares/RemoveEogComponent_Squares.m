function EEGout = RemoveEogComponent_Squares(EEG,sq,samples,component)

% Created 7/2/12 by DJ.

% Get pixel coordinates of center of screen
Constants = GetSquaresConstants;
X_CENTER = Constants.SQUARE_X(3);
Y_CENTER = Constants.SQUARE_Y(3);

% Interpolate eye position at time of each EEG sample
t = (1:size(samples,1)) + sq.recording_start_time - 1;
tEeg = EyelinkToEegTimes(t,sq);
eyedata = interp1(tEeg,samples(:,1),1:EEG.pnts,'linear',NaN) - X_CENTER;

% Fill in blink times (but leave times outside recording as NaNs)
isNanTime = isnan(eyedata);
eyedata(isNanTime) = interp1(find(~isNanTime),eyedata(~isNanTime),find(isNanTime),'linear',NaN);

% Predict EOG contribution to signal
EOG_prediction = component'*eyedata;
% Turn the nans into zeros
EOG_prediction(isnan(EOG_prediction)) = 0;

% Prepare output
EEGout = EEG;
EEGout.data = EEG.data-EOG_prediction;


