function WriteTspTextFile(subject,sessions,levelname,TspTour)

% Created 12/13/12 by DJ.

% Get level info
levelinfo = load([levelname(1:end-4) '.mat']); % load variables campoints and sessionOffset

% main loop
[objects, ~, ~, objtimes, objisessions] = GetObjectList(subject,sessions); 
objoffset = nan(length(objects),2);
for i=1:numel(sessions)
    thisSessionOffset = levelinfo.sessionOffset * (i-1);
    objoffset(objisessions==i,:) = repmat(thisSessionOffset,sum(objisessions==i),1);
end
% Make string for whether each object was seen
objseen(~isnan(objtimes)) = {'Seen'};
objseen(isnan(objtimes)) = {'Unseen'};

% Create text file
foonew = fopen(sprintf('3DS-%d-tspout.txt',subject),'w');

fprintf(foonew,'----- SESSION PARAMETERS -----\n');
% Copy session parameters from an input text file
foo = fopen(sprintf('3DS-%d-%d.asc',subject,sessions(1)),'r');
thisline = fgetl(foo);
while isempty(strfind(thisline,'SESSION PARAMETERS'))
    thisline = fgetl(foo);
end
thisline = fgetl(foo);
while isempty(strfind(thisline,'END SESSION PARAMETERS'))    
    fprintf(foonew,'%s\n',thisline);
    thisline = fgetl(foo);
end
fclose(foo);

fprintf(foonew,'----- END SESSION PARAMETERS -----\n');

% Write object info
fprintf(foonew,'----- OBJECT INFO -----\n');
for i=1:length(objects)
    fprintf(foonew,'Object # %d %s %s %s (%g %g %g) (%g %g %g %g) %s\n',...
        i,objects(i).name,objects(i).tag,objects(i).type,...
        objects(i).createposition(1)+objoffset(i,1),objects(i).createposition(2),objects(i).createposition(3)+objoffset(i,2),...
        objects(i).createrotation(1),objects(i).createrotation(2),objects(i).createrotation(3),objects(i).createrotation(4), objseen{i});
%     EXAMPLE: fprintf(foonew,'Object # 1 schooner-24 DistractorObject Stationary (94.25, 0.75, 141) (0, 1, 0, -4.371139E-08)');
end
fprintf(foonew,'----- END OBJECT INFO -----\n');

% Write TSP route info
fprintf(foonew,'----- ROUTE -----\n');
for i=1:size(TspTour,1)
    fprintf(foonew,'Point #%d (%d,%d)\n',i,TspTour(i,1),TspTour(i,2));
end
fprintf(foonew,'----- END ROUTE -----\n');


fclose(foonew);

