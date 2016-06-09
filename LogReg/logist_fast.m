function [ws_loo,ys_loo,Az_loo,Az_perm_dist] = logist_fast(X,y,l2_lambda,mode,varargin)
% X is the feature data (matrix size DxN -- D = # features, N = # trials)
% y are the binary truth labels (vector size Nx1)
% l2_lambda is the l2 regularization value (e.g., 1e-6)
% mode can be either 'loo' or 'permutation'

X = double(X);
y = double(y);
l2_lambda = double(l2_lambda);
conv_tol = 1e-11;

looMode = 'average';

% Add all one's feature -- to estimate the bias
[D,N] = size(X);
X = [X;ones(1,N)];
D = D+1;
debug = 0;

chunkSize = 5000; % Can increase this number if it doesn't crash!
numBatches = 1;
if strcmp(mode,'loo')
    if isempty(varargin)
        cv = [];
        cv.numFolds = N;
        cv.outTrials = cell(1,cv.numFolds);
        for p=1:cv.numFolds
            cv.outTrials{p} = p;
        end
    else
        cv = varargin{1};
    end
    
    yperms = assertVec(y,'col');    
elseif strcmp(mode,'loo_and_bootstrap')
    cv = varargin{1};
    yperms = [assertVec(y,'col'),varargin{2}];
else    
    error('Unknown mode');
end
    
numperms = size(yperms,2);
% Make chunkSize the next closest multiple of cv.numFolds
chunkSize = cv.numFolds*round(chunkSize/cv.numFolds);
P = cv.numFolds*numperms;
if P > chunkSize
    numBatches = ceil(P/chunkSize);
    P = chunkSize;
    disp(['There are too many permutations to run in one batch.  This analysis will be run in ',num2str(numBatches), ' batches.']);pause(1e-9);
end
    
[y,zeroOutInds,perms] = getCurrentYAndZeroInds(1,yperms,cv,chunkSize);

ws_loo = zeros(D,cv.numFolds);
if strcmp(looMode,'average')
    ys_loo = zeros(cv.numFolds,1);
else
    ys_loo = zeros(N,1);
end
Az_perm_dist = zeros(numperms-1,1);

% Set the regularization matrix L
% Set the last diagonal element to zero to not penalize the bias
L = diag([ones(D-1,1);0]);

[Q,Z] = qr(X,0);
if size(Q,2) < D
    % Then we're better off using a low-rank representation for X
      
    % Now we need to find a basis for components y for which L*y is in the
    % range of U.  Also the component y must be orthogonal to U (U'*y=0)
    [Qperp,~] = qr(L'*Q,0);
    Qperp = (eye(size(Q,1))-Q*Q')*Qperp;
    [Qtmp,~] = qr([Q,Qperp],0);
    Q = Qtmp(:,1:size(Q,2));
    Qperp = Qtmp(:,(size(Q,2)+1):end);
    
    B = Q'*L*Qperp;
    D = Qperp'*L*Qperp;
    DinvBT = D\(B'); % DinvBT = D^(-1)*B';
    
    additiveMatrixFactor = 2*l2_lambda*(Q'*L*Q - B*DinvBT);
    
    D = size(Q,2);
    finalInvertMat = Q - Qperp*DinvBT;
    
    error('Didn''t reset Z here!');
else
    additiveMatrixFactor = 2*l2_lambda*L;
    finalInvertMat = eye(D);
    Z = X;
end

XT = X';

innerLoopMode = 'singleMat';'useG';
if strcmp(innerLoopMode,'useG')
    G = cell(1,P);
    for p=1:P
        G{p} = zeros(D,D);
    end
end

ws = zeros(D,P);
% Pre-allocate some of the big matrices
ws_last = zeros(D,P);
wsprev = ws;
onesMat = ones(N,P);
RTREP = zeros(N,P);
RTMAT = sparse(N,N);
RT_minus_R = zeros(N,P);
yvals = zeros(N,P);
mus = zeros(N,P);
R = zeros(N,P);
Minv_ls = zeros(D,P);

% OK now let's get started with the iterations
iter = 0;
locsNotFullyConverged = 1:P;
batchNum = 1;
while 1
    iter = iter + 1;
    disp(['Working on batch #',num2str(batchNum),' out of ',num2str(numBatches),'; Iteration #',num2str(iter),'...']);pause(1e-9);

    ws_last(:,:) = ws;

    yvals(:,locsNotFullyConverged) = XT*ws(:,locsNotFullyConverged);
	mus(:,locsNotFullyConverged) = 1./(1+exp(-yvals(:,locsNotFullyConverged)));

    R(:,locsNotFullyConverged) = mus(:,locsNotFullyConverged).*(1-mus(:,locsNotFullyConverged));
    R(zeroOutInds) = 0;
                
	RT = max(R(:,locsNotFullyConverged),[],2);%%median(R(:,locsNotFullyConverged),2);
	if max(RT) == 0
	    break;
    end
        
	neg_ls_without_X(:,locsNotFullyConverged) = (R(:,locsNotFullyConverged).*yvals(:,locsNotFullyConverged) + y(:,locsNotFullyConverged) - mus(:,locsNotFullyConverged));
    neg_ls_without_X(zeroOutInds) = 0;
%    neg_ls = X*neg_ls_without_X;

    RTMAT = spdiags(RT,0,RTMAT);
    M = X*(RTMAT*XT) + additiveMatrixFactor;
    MinvX = M\X;
	if strncmp(lastwarn,'Matrix is close to singular or badly scaled.',44)
	    warning('Bad conditioning. Suggest reducing subspace.')
	    singularwarning=1;
	    break;
    end
                
	% Iteration is the following:
	% ws = Minv*N*wsprev + Minv*neg_ls
	Minv_ls(:,locsNotFullyConverged) = MinvX*neg_ls_without_X(:,locsNotFullyConverged);
    RTREP(:,locsNotFullyConverged) = RTMAT*onesMat(:,locsNotFullyConverged);
%    RTREP = repmat(RT,1,length(locsNotFullyConverged));
	RT_minus_R(:,locsNotFullyConverged) = RTREP(:,locsNotFullyConverged) - R(:,locsNotFullyConverged);        
        
    
    performLowRankUpdates = 0;
    ratio_thresh = 0.5;
    maxRank = 1;%round(0.05*D);

    hasLowRankCorrections = zeros(1,P);
    lowRankLocs = cell(1,P);
    lowRankCorrectionMatsLeft = cell(1,P);
    lowRankCorrectionMatsRight = cell(1,P);
    lowRankCorrectionSizes = zeros(1,P);
    if performLowRankUpdates == 1
        RT_ratio = zeros(N,P);    
        RT_ratio(:,locsNotFullyConverged) = abs(RT_minus_R(:,locsNotFullyConverged))./RTREP(:,locsNotFullyConverged);
        RT_ratio(zeroOutInds) = 0;disp('Zeroing out RT_ratio(zeroOutInds) for now.  May not want to do this in the future to speed up the algorithm.');pause(1e-9);
        colinds = max(RT_ratio(:,locsNotFullyConverged)) > ratio_thresh;
        lowRankColIndsUnique = locsNotFullyConverged(colinds);    
        avgLowRankNum = 0;
        if ~isempty(lowRankColIndsUnique)
            XTMinvX = XT*MinvX;
            for p=lowRankColIndsUnique
                lowRankLocs{p} = find(RT_ratio(:,p) > ratio_thresh);
                if isempty(lowRankLocs{p}); 
                    continue; 
                end;
    
                if length(lowRankLocs{p}) > maxRank
                    % Then make sure that we get everywhere that RT_ratio == 1
                    numToKeep = max(maxRank,length(find(RT_ratio(lowRankLocs{p},p) == 1)));
                    [~,inds] = sort(RT_ratio(lowRankLocs{p},p),'descend');
                    lowRankLocs{p} = lowRankLocs{p}(inds(1:numToKeep));
                end

                lowRankCorrectionSizes(p) = length(lowRankLocs{p});
                avgLowRankNum = avgLowRankNum + length(lowRankCorrectionSizes(p));

                hasLowRankCorrections(p) = 1;
                lowRankCorrectionMatsLeft{p} = MinvX(:,lowRankLocs{p});
                lowRankCorrectionMatsRight{p} = (diag(1./RT_minus_R(lowRankLocs{p},p)) - XTMinvX(lowRankLocs{p},lowRankLocs{p}))\XT(lowRankLocs{p},:);%XT(lowRankLocs{p},:)*MinvX(:,lowRankLocs{p})

                RT_minus_R(lowRankLocs{p},p) = 0;
                RT_ratio(lowRankLocs{p},p) = 0;        
            end
        end
    else
        disp('Not performing low rank updates.  May want to do this in the future to speed up the algorithm.');pause(1e-9);
    end
    hasLowRankCorrectionsCell = num2cell(hasLowRankCorrections);
%    avgLowRankNum/length(locsNotFullyConverged)
    
    if strcmp(innerLoopMode,'useG')
        for p=locsNotFullyConverged
            G{p}(:,:) = MinvX*(repmat(RT_minus_R(:,p),1,D).*XT);
        end        
        Minvbc = mat2cell(Minv_ls,D,ones(1,P));
    end
        
	locsNotConverged = locsNotFullyConverged;
	ninnerloops = 0;
	while 1
        ninnerloops = ninnerloops + 1;
        %fprintf(1,['\tPerforming inner loop #',num2str(ninnerloops),'...\n']);pause(1e-9);
        wsprev(:,:) = ws;

        if strcmp(innerLoopMode,'useG')
            wsc = mat2cell(ws(:,locsNotConverged),D,ones(1,length(locsNotConverged)));
            ws(:,locsNotConverged) = cell2mat(cellfun(@iterateWithLowRankCorrection,wsc,G(locsNotConverged),Minvbc(locsNotConverged),lowRankCorrectionMatsLeft(locsNotConverged),lowRankCorrectionMatsRight(locsNotConverged),hasLowRankCorrections(locsNotConverged),'UniformOutput',false));
        else
            %XTwsprev = XT*wsprev(:,locsNotConverged);
            ws(:,locsNotConverged) = MinvX*(RT_minus_R(:,locsNotConverged).*(XT*wsprev(:,locsNotConverged))) + Minv_ls(:,locsNotConverged);

            % Now do the low-rank corrections
            lowRankLocsCurr = locsNotConverged(find(hasLowRankCorrections(locsNotConverged)));
            if ~isempty(lowRankLocsCurr)
                wsc = mat2cell(ws(:,lowRankLocsCurr),D,ones(1,length(lowRankLocsCurr)));
                ws(:,lowRankLocsCurr) = cell2mat(cellfun(@lowRankCorrection,wsc,lowRankCorrectionMatsLeft(lowRankLocsCurr),lowRankCorrectionMatsRight(lowRankLocsCurr),hasLowRankCorrectionsCell(lowRankLocsCurr),'UniformOutput',false));
            end
        end

        tol_measure = (1/D)*sum((ws(:,locsNotConverged)-wsprev(:,locsNotConverged)).^2)./sum(ws(:,locsNotConverged).^2);
        locsNew = find(tol_measure > conv_tol);
        if isempty(locsNew)
            break;
        end
        locsNotConverged = locsNotConverged(locsNew);

        if debug==1
            figure(102);
            plot((sum((ws-wsprev).^2)./sum(wsprev.^2)));
            pause(0.1);
        end
    end
    fprintf(1,'\tNumber of inner loops: %d\n',ninnerloops);pause(1e-9);
        
    locsNotFullyConverged = find((1/D)*sum((ws_last-ws).^2)./sum(ws_last.^2) > conv_tol);
    
    if isempty(locsNotFullyConverged)
        % Save out the converged results
        eInd = 0;
        while 1
            sInd = eInd + 1;
            if sInd > length(perms)
                break
            end
            eInd = sInd + cv.numFolds - 1;
            currInds = sInd:eInd;
            if max(abs(perms(currInds)-mean(perms(currInds))))
                error('Something''s very wrong here!');
            end
            permInd = perms(sInd);

            if strcmp(looMode,'average')
                ys_loo_perm = zeros(cv.numFolds,1);
                truth = zeros(cv.numFolds,1);
            else
                ys_loo_perm = zeros(N,1);
                truth = zeros(N,1);
            end
                
            for foldNum = 1:cv.numFolds
                if strcmp(looMode,'average')
                    ys_loo_perm(foldNum) = mean(XT(cv.outTrials{foldNum},:)*ws(:,currInds(foldNum)));
                    truth(foldNum) = yperms(cv.outTrials{foldNum}(1),permInd);
                else
                    ys_loo_perm(cv.outTrials{foldNum}) = XT(cv.outTrials{foldNum},:)*ws(:,currInds(foldNum));
                    truth(cv.outTrials{foldNum}) = yperms(cv.outTrials{foldNum},permInd);
                end
            end
            % Now compute the roc value and store it
            Az_val = rocarea(ys_loo_perm,truth);
            
            if permInd == 1
                % Then this is the LOO run -- save out the ws
                ws_loo(:,:) = ws(:,currInds);
                ys_loo = ys_loo_perm;
                Az_loo = Az_val;
            else
                Az_perm_dist(permInd-1) = Az_val;
            end
        end
        
        if batchNum < numBatches        
            % We need to re-fill
            iter = 0;
            batchNum = batchNum + 1;
            ws = zeros(D,P);
            [y,zeroOutInds,perms,locsNotFullyConverged] = getCurrentYAndZeroInds(batchNum,yperms,cv,chunkSize);
        else
            break;
        end
    end
end
%toc
ws = finalInvertMat*ws;


%if strcmp(mode,'loo')
%    yv = sum(X.*ws)';
%    testerr = 0.5*sum(abs((2*(y(:,1)-0.5)) - sign(yv)));
%    disp(['Number of test errors:  ',num2str(testerr)]);
%end

%disp(['DONE!; ', num2str(iter)]);
end


function w = lowRankCorrection(w,lowRankCorrectionMatLeft,lowRankCorrectionMatRight,hasLowRankCorrection)
    if ~hasLowRankCorrection; return; end;
    
    w = w + lowRankCorrectionMatLeft*(lowRankCorrectionMatRight*w);
end

function e = lowRankCorrection2(w,lowRankCorrectionMatLeft,lowRankCorrectionMatRight,hasLowRankCorrection)
    e = zeros(size(w),class(w));
    if ~hasLowRankCorrection; return; end;
    
    e = lowRankCorrectionMatLeft*(lowRankCorrectionMatRight*w);
end

function w = iterateWithLowRankCorrection(w,G,Minvb,lowRankCorrectionMatLeft,lowRankCorrectionMatRight,hasLowRankCorrection)
    w = G*w + Minvb;
    if hasLowRankCorrection
        w = w + lowRankCorrectionMatLeft*(lowRankCorrectionMatRight*w);
    end
end



function [y,zeroOutInds,perms,locsNotFullyConverged] = getCurrentYAndZeroInds(batchNum,yperms,cv,chunkSize)
    [N,numperms] = size(yperms);    
    P = min(chunkSize,cv.numFolds*numperms);
    
    numpermsPerBatch = P/cv.numFolds;
    if abs(numpermsPerBatch - round(numpermsPerBatch)) > 0
        error('Something''s not right');
    end
    % Now determine the set of permutations to run based on batchNum
    permsStart = (batchNum-1)*numpermsPerBatch + 1;
    permsEnd = min(numperms, permsStart + numpermsPerBatch - 1);
    permsCurr = permsStart:permsEnd;
    numpermsCurr = length(permsCurr);
    
    perms = reshape(repmat(assertVec(permsCurr,'row'),cv.numFolds,1),numpermsCurr*cv.numFolds,1);
    
    y = zeros(N,P);
    eInd = 0;
    for n=permsCurr
        sInd = eInd + 1;
        eInd = sInd + cv.numFolds - 1;
        y(:,sInd:eInd) = repmat(yperms(:,n),1,cv.numFolds);
    end
    
    locsNotFullyConverged = 1:eInd;
    
    totalOutTrials = 0;
    for p=1:cv.numFolds
        totalOutTrials = totalOutTrials + length(cv.outTrials{p});
        cv.outTrials{p} = assertVec(cv.outTrials{p},'row');
    end
    
    zeroOutInds = zeros(1,totalOutTrials*numpermsCurr,'uint32');
    eInd = 0;
    % Compute the indices to be zeroed out in the NxP gradient matrix
    for foldNum=1:cv.numFolds
        for n=1:numpermsCurr
            sInd = eInd + 1;
            eInd = sInd + length(cv.outTrials{foldNum}) - 1;
            zeroOutInds(sInd:eInd) = cv.outTrials{foldNum} + N*(foldNum+(n-1)*cv.numFolds - 1);%sub2ind([N,P],assertVec(cv.outTrials{foldNum},'row'),(foldNum+(n-1)*cv.numFolds)*ones(1,length(cv.outTrials{foldNum})));
        end
    end
end


function doSomePermutationStuff
            ratio_thresh = 0.99;
            max_inversions = 10;%0.1*P;
            
            % We need to find extreme RT_minus_R/RT values (close to 1)
            % that may slow convergence
            QT = RT_minus_R(:,locsNotFullyConverged)./repmat(RT,1,length(locsNotFullyConverged));
            [nctrialinds,ncperminds] = find(QT > ratio_thresh);
            ncperminds = locsNotFullyConverged(ncperminds);
            uncperminds = unique(ncperminds);
            [~,permindordering] = sort(min(R(:,uncperminds)),'ascend');
            uncperminds = [];%uncperminds(permindordering);
            
            numclusts = 1;            
            cluster_perminds = cell(1,0);
            cluster_perminds{1} = setdiff(locsNotFullyConverged,uncperminds);
            RT_minus_Rs = cell(1,0);
            RT_minus_Rs{1} = repmat(RT,1,length(cluster_perminds{1})) - R(:,cluster_perminds{1});
            Minvs = cell(1,0);
            Minvs{1} = Minv;
            MinvXs = cell(1,0);
            MinvXs{1} = MinvX;
            Minv_lsc = cell(1,0);
            Minv_lsc{1} = Minv_ls(:,cluster_perminds{1});
            if length(uncperminds)
                while length(uncperminds)
                    uncpermind = uncperminds(1);
                    locs = nctrialinds(find(ncperminds==uncpermind));
                    RTC = RT;
                    RTC(locs) = R(locs,uncpermind)/(1-(ratio_thresh-eps));
                    
%                    locsprev = [];
%                    while 1
                        locs = uncperminds(find(min(repmat(RTC,1,length(uncperminds))-R(:,uncperminds))>=0));
                        RTC = max(R(:,locs),[],2);                        
%                        locs = intersect(locs,find(min(R(:,uncperminds)/(1-(ratio_thresh-eps))-repmat(RTC,1,length(uncperminds)))>=0));
                        
%                        if (length(locs) == length(locsprev)) && ~length(setdiff(locs,locsprev))
%                            break;
%                        end
%                        locsprev = locs;
%                    end
%                    if max(max((repmat(RTC,1,length(locs))-R(:,locs))./repmat(RTC,1,length(locs)))) > ratio_thresh
%                        error('Something wrong here!');
%                    end
                    numclusts = numclusts + 1;
                    cluster_perminds{numclusts} = locs;
                    
                    ismem = ismember(uncperminds,cluster_perminds{numclusts});
                    uncperminds = uncperminds(~ismem);
%                    disp([num2str(numclusts),'; ',num2str(length(uncperminds)),'; ',num2str(max(max((repmat(RTC,1,length(cluster_perminds{numclusts}))-R(:,cluster_perminds{numclusts}))./repmat(RTC,1,length(cluster_perminds{numclusts})))))]);pause(1e-9);
                    
                    RT_minus_Rs{numclusts} = repmat(RTC,1,length(cluster_perminds{numclusts}))-R(:,cluster_perminds{numclusts});
                    Minvs{numclusts} = (X*diag(RTC)*XT + additiveMatrixFactor)^(-1);
                    MinvXs{numclusts} = Minvs{numclusts}*X;
                    Minv_lsc{numclusts} = MinvXs{numclusts}*neg_ls_without_X(:,cluster_perminds{numclusts});
                    
                    if numclusts == (max_inversions-1)
                        numclusts = numclusts + 1;
                        cluster_perminds{numclusts} = uncperminds;
                        RTC = max(R(:,uncperminds),[],2);
                        RT_minus_Rs{numclusts} = repmat(RTC,1,length(cluster_perminds{numclusts}))-R(:,cluster_perminds{numclusts});
                        Minvs{numclusts} = (X*diag(RTC)*XT + additiveMatrixFactor)^(-1);
                        MinvXs{numclusts} = Minvs{numclusts}*X;
                        Minv_lsc{numclusts} = MinvXs{numclusts}*neg_ls_without_X(:,cluster_perminds{numclusts});
                        break;
                    end
                end
            else
                % Do nothing then...
            end
%            numclusts
end


% 	if strcmp(mode,'permutation')
% %                ws_norc = ws;
%                 turnNum = 0;
%                 numTurns = 1;
%                 cachesize = 3;
%                 ws_cache = zeros(D,size(ws,2)*cachesize,class(ws));
%                 while 1
%                     turnNum = turnNum + 1;
%                     if turnNum == numTurns
%                         conv_tol = 1e-12;
%                     else
%                         conv_tol = 1e-6;
%                     end
%                     
%                     locsNotConverged = locsNotFullyConverged;
%                     iterc = 0;
%                     while 1
%                         iterc = iterc + 1;
%                         wsprev = ws; 
%                         
% %                        XTwsprev = XT*wsprev(:,locsNotConverged);                    
%                         ws(:,locsNotConverged) = MinvX*(RT_minus_R(:,locsNotConverged).*(XT*wsprev(:,locsNotConverged))) + Minv_ls(:,locsNotConverged);
%                         if turnNum == numTurns
%                             % Now do the low-rank corrections
%                             wsc = mat2cell(ws(:,locsNotConverged),D,ones(1,length(locsNotConverged)));
%                             ws(:,locsNotConverged) = cell2mat(cellfun(@lowRankCorrection,wsc,lowRankCorrectionMatsLeft(locsNotConverged),lowRankCorrectionMatsRight(locsNotConverged),hasLowRankCorrections(locsNotConverged),'UniformOutput',false));
%                         else
%                             for c=min(iterc,cachesize):-1:1
%                                 inds = c:cachesize:size(ws_cache,2);
%                                 inds = inds(locsNotConverged);
%                                 if c > 1
%                                     indsprev = (c-1):cachesize:size(ws_cache,2);
%                                     indsprev = indsprev(locsNotConverged);
%                                     ws_cache(:,inds) = ws_cache(:,indsprev);
%                                 else
%                                     ws_cache(:,inds) = ws(:,locsNotConverged);
%                                 end
%                             end
%                         end
%                         
%                         tol_measure = (1/N)*sum((ws(:,locsNotConverged)-wsprev(:,locsNotConverged)).^2)./sum(ws(:,locsNotConverged).^2);
%                         %tol_measure = max(abs(ws(:,locsNotConverged)-wsprev(:,locsNotConverged))./abs(ws(:,locsNotConverged)));
%                         locsNew = find(tol_measure > conv_tol);
%                         if isempty(locsNew)
%                             break;
%                         end
%                         locsNotConverged = locsNotConverged(locsNew);
%                 
%                         if debug==1
%                             h=figure(102);
%                             plot(tol_measure);
%                         end
%                     end
%                     iterc
%                     
% %                    wstrue = zeros(size(ws));
% %                    for p=1:P
% %                        wstrue(:,p) = (X*diag(R(:,p))*XT+2*l2_lambda*L)\neg_ls(:,p);
% %                    end
% %                    figure(107);plot((1/N)*sum((ws-wstrue).^2)./sum(wstrue.^2));
% %                    disp('DONE!!!');
% %                    pause;
%                     
%                     if turnNum == numTurns
%                         break;
%                     end
%                     
%                     e = zeros(D,P);
%                     wscc = mat2cell(ws_cache,D,cachesize*ones(1,P));
%                     es = zeros(D,P*cachesize);
%                     inds = reshape(1:P*cachesize,cachesize,P);
%                     inds = reshape(inds(:,locsNotFullyConverged),1,cachesize*length(locsNotFullyConverged));
%                     es(:,inds) = cell2mat(cellfun(@lowRankCorrection2,wscc(locsNotFullyConverged),lowRankCorrectionMatsLeft(locsNotFullyConverged),lowRankCorrectionMatsRight(locsNotFullyConverged),hasLowRankCorrections(locsNotFullyConverged),'UniformOutput',false));
%                     for cacheind = min(iterc,cachesize):-1:1
%                         inds = cacheind:cachesize:size(es,2);
%                         inds = inds(locsNotFullyConverged);
%                         e(:,locsNotFullyConverged) = e(:,locsNotFullyConverged) + es(:,inds);
%                         if cacheind > 1
%                             e(:,locsNotFullyConverged) = MinvX*(RT_minus_R(:,locsNotFullyConverged).*(XT*e(:,locsNotFullyConverged)));
%                             e(:,locsNotFullyConverged) = cell2mat(cellfun(@lowRankCorrection,mat2cell(e(:,locsNotFullyConverged),D,ones(1,length(locsNotFullyConverged))),lowRankCorrectionMatsLeft(locsNotFullyConverged),lowRankCorrectionMatsRight(locsNotFullyConverged),hasLowRankCorrections(locsNotFullyConverged),'UniformOutput',false));
%                         end
%                     end
%                     ws = ws + e;
% %                        figure(108);plot((1/N)*sum((ws-wstrue).^2)./sum(wstrue.^2));
% %                        disp(num2str(cacheind));
% %                        pause;
%                     
%                 end
% %            end
% %    end

