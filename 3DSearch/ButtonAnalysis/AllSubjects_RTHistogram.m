% This script can be used to get the RT histograms for each subject (and
% combined across subjects), and the diffusion model parameters that best
% fit these histograms.
%
% AllSubjects_RTHistogram.m
%
% It's a bit messy because it was designed in a hurry. 
% TO DO: Clean it up!
%
% Created 3/24/11 by DJ.

subjects = [8,10,11,12];
% Get RTs
[~,~,RT{1}] = CompiledRTHistogram(8,[2:5,7:9 11],'objects');

[~,~,RT{2}] = CompiledRTHistogram(10,[2:3, 5:10],'objects');

[~,~,RT{3}] = CompiledRTHistogram(11,[2:5,7:10],'objects');

[~,~,RT{4}] = CompiledRTHistogram(12,3:10,'objects');
AllRT = [RT{:}];

nSubj = numel(RT);

% Remove outliers
for i=1:nSubj
    RT{i}(RT{i}>mean(RT{i})+3*std(RT{i})) = [];
end

mu = mean(AllRT);
sigma = std(AllRT);
AllRT(AllRT>mu+3*sigma) = [];

% Plot
clf
t = 300:50:1000;
for i=1:nSubj
    subplot(4,2,i)
    n{i} = hist(RT{i},t);
    bar(t,n{i});
end

subplot(2,1,2)
nAll = hist(AllRT,t);
bar(t,nAll)

%% Make diffusion model pdf with hand-picked parameters
L = 600;
r = 1.1;
d = 4;

% pdf = L ./ sqrt(2*pi*d^2 * t.^3) .* exp(-(L-r*t).^2 ./ (2*d^2*t));
% pdf = pdf/sum(pdf)*sum(nAll);
% 
% % Plot over desired axis
% subplot(2,1,2)
% hold on;
% plot(t,pdf,'g');
% 
% for i=1:nSubj
%     subplot(4,2,i); hold on;
%     plot(t,pdf/sum(pdf)*numel(RT{2}),'g');
% end

%% Get min-mean-squared-error pdf for each subject

% Get ideal plot for each individual subject
for i=1:nSubj
    % use fminsearch
    [x, mse] = fminsearch(@(x) pdf_mse(x,t,n{i}),[L r d]);
    % parse results
    L = x(1);
    r = x(2);
    d = x(3);
    
    % compute pdf
    pdf = L ./ sqrt(2*pi*d^2 * t.^3) .* exp(-(L-r*t).^2 ./ (2*d^2*t));
    pdf = pdf/sum(pdf)*sum(n{i});

    % Plot over desired axis
    subplot(4,2,i)
    hold on;
    plot(t,pdf,'r');
    % Display results
    fprintf('Subject #%d: L=%.0f, r=%.3g, d=%.3g, mse=%.3g\n',subjects(i),L,r,d,mse);
    title(sprintf('Subject #%d: L=%.0f, r=%.3g, d=%.3g, mse=%.3g\n',subjects(i),L,r,d,mse));
end


% ...And do the same for the aggregated data
% use fminsearch
x = fminsearch(@(x) pdf_mse(x,t,nAll),[L r d]);
% parse results
L = x(1);
r = x(2);
d = x(3);
% compute pdf
pdf = L ./ sqrt(2*pi*d^2 * t.^3) .* exp(-(L-r*t).^2 ./ (2*d^2*t));
pdf = pdf/sum(pdf)*sum(nAll);

% Plot over desired axis
subplot(2,1,2)
hold on;
plot(t,pdf,'r');
% Display results
fprintf('All Subjects: L=%.0f, r=%.3g, d=%.3g, mse=%.3g\n',L,r,d,mse);
title(sprintf('All Subjects: L=%.0f, r=%.3g, d=%.3g, mse=%.3g\n',L,r,d,mse));

