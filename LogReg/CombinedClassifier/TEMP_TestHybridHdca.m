N = 100;
T = 100;
D = 50;
truth = round(rand(1,N)); % set random trials to be true

twl = 10;
two = 1:10:100;

data = rand(D,T,N);
data(1:20,31:50,truth==1) = data(1:20,31:50,truth==1) + .5;
data = data - repmat(mean(data,3),[1,1,N]);

level2data = rand(N,2);
level2data(truth==1,1) = level2data(truth==1,1) + 1;

% [y, w, v, fwdModel, y_level1] = TrainHybridHdcaClassifier(data, truth, twl, two, level2data);

% Compare CV classifiers
% [yCV, wCV, vCV, fwdModelCV, y_level1CV] = RunHdcaClassifier(data, truth, twl, two, '10fold');
[yCV2, wCV2, vCV2, fwdModelCV2, y_level1CV2] = RunHybridHdcaClassifier(data, truth, twl, two, '10fold',level2data);


%% PLOT
figure
% subplot(211);
% imagesc(w);
% subplot(212);
% plot(v);
subplot(221);
imagesc(mean(mean(wCV,3),4));
subplot(222);
imagesc(mean(mean(wCV2,3),4));
subplot(223);
plot(mean(vCV,3));
subplot(224);
plot(mean(vCV2,3));

%% TEST

testdata = data;
testtruth = truth;
% testtruth(1:20) = 0;

yTest = ApplyHybridHdcaClassifier(testdata,twl,two,level2data,w,v);
AzTest = rocarea(yTest,testtruth);
fprintf('AzTest = %.2f\n',AzTest);