%%% 08/24/2007   Jianing Shi
%%% Analyze classification results
%%% Classification option
%%% option 1 : naive linear classification error, percent of correctness
%%% option 2 : Az value through ROC analysis

function res = ClassificationResult(X, b, w, v, opt)

%%% Naive
if opt == 1

proj = X*w + v;
pred = sign(proj);
totalN = length(b);
corrN = sum(pred == b);
res = corrN/totalN;

%%% ROC
elseif opt == 2

truth = (b == 1);
proj = X*w + v;
bp = bernoull(1,proj);
res = rocarea(bp, truth);

%%% ROC using auc
elseif opt == 3
    
truth = (b == 1);
proj = X*w + v;
bp = bernoull(1,proj);
res = auc(truth, bp);

%%% Plot error
elseif opt == 4
    
proj = X*w + v;
pred = sign(proj);
totalN = length(b);
errorN = sum(pred ~= b);
res = errorN/totalN;   
    
end