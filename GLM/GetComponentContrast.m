function GetComponentContrast(rules,experiments,old_suffix,suffix_in,suffix_out,iEvents_cell,lambda_best,betas,components)

% GetComponentContrast(rules,experiments,old_suffix,suffix_in,suffix_out,iEvents_cell,lambda_best,betas,components)
%
% INPUTS:
% -components is an ExD matrix.
%
% Created 9/19/14 by DJ.
% Updated 1/15/15 by DJ -changed to iEvents_cell


if ~exist('components','var') || isempty(components)
    components = eye(size(betas{1},1));    
end

multcorrect = 'none';

for iRule = 1:numel(rules)
    rule = rules{iRule};
    for iExp=1:numel(experiments)
        experiment = experiments{iExp};    
        [subjects,basedir,folders] = GetSquaresSubjects(experiment);
        cd(basedir);
        fprintf('=== %s %s, %s: Setting up...\n',datestr(now,16),rule,experiment);
        % set up contrast
        R = load(sprintf('%s-%d-%s',experiment,subjects(1),old_suffix));
        iEvents = iEvents_cell{iExp};
        
        p = length(R.tResponse{end})*numel(iEvents);
        event_list = R.regressor_events{end}(iEvents);
        tResponse = R.tResponse{end};
        chanlocs = R.EEG.chanlocs;
        
        % Get contrast matrix
        [contrasts0,legendstr,sequence,tContrast] = GetSequenceContrast(experiment,rule,event_list,tResponse);


        %--- Get Single-Subject Contrasts
        lambda = lambda_best(iExp);
        % set up contrast functions
        D = size(components,1);%numel(chanlocs);
        N = numel(subjects);    
        T = length(tResponse);
        M = size(contrasts0,2)/T; % number of contrasts
        [contrastFns, contrastVar, contrastZ] = deal(nan(D,T,M,N));

        % find contrast functions
        for iSubj=1:N
            fprintf('---subject %d/%d...\n',iSubj,N); 
            file_in = sprintf('%s-%d-%s',experiment,subjects(iSubj),suffix_in);
            R = load(file_in);
            R.Y = R.Y*components';
            
            beta = components*betas{iExp}(:,:,iSubj);
            [sigmasq_all] = GetGlmMse(R.X,R.Y,beta);
            
            
            XTX = full(R.X'*R.X);
%             XTY = full(R.X')*R.Y;
            n = size(R.X,1);
            Xmean = R.Xmean;
            Xrss = R.Xrss;
            clear R;

            % Calculate and plot
        %     [contrastZ(:,:,:,i),~,~,contrastFns(:,:,:,i),contrastVar(:,:,:,i)] = GetGlmZscore(R(i).EEG,R(i),multcorrect,contrasts,doPlot);
            fprintf('Calculating dof... ');    
            E = eig(XTX);
            dof = n - sum(E./(E+lambda));  % df(lambda) = sum(i=1:p)(d_i^2/(d_i^2+lambda))
            fprintf('dof = %g. \nCalcluating inverse...\n',dof);
            tic;
            X_inv = (XTX + lambda*eye(p))^(-1);


            contrasts = contrasts0;
            for i=1:p
                contrasts(i,:) = (contrasts(i,:)-Xmean(i))/Xrss(i);
            end


            


            fprintf('Getting correction factors...')
            tic;
            % Get correction factors
            correction_fact = nan(D,M*T);
            for i=1:(M*T)
                correction_fact(:,i) = contrasts(:,i)'*X_inv*contrasts(:,i); % correction factor to account for linear dependence between contrasts
            end
            t=toc;
            fprintf('Done! Took %.1f seconds.\n',t);

            fprintf('Calculating contrasts...')
            tic;
        %     contrastFns = zeros(D,T,M); % like responseFns
        %     contrastVar = zeros(D,T,M);
            for j=1:D
                for i=1:size(contrasts,2)    
                    % Get inidces
                    iCF = ceil(i/T); % which contrast fn?
                    iT = i-T*(iCF-1); % which time point in that CF?

            %         contrastVar(:,iT,iCF) = sigmasq'*correction_fact(:,i)/dof; % variance of contrast function
                    contrastVar(j,iT,iCF,iSubj) = sigmasq_all*correction_fact(j,i)/dof; % variance of contrast function

                    % Get response fn of contrast    
            %         contrastFns(:,iT,iCF) = beta*contrasts(:,i);
                    contrastFns(j,iT,iCF,iSubj) = beta(j,:)*contrasts(:,i); % convert response functions to contrast functions
                end    
            end
            t=toc;
            fprintf('Done! Took %.1f seconds.\n',t);
        end
        disp('Getting statistics...')
        tic
        % take sqrt to get std error (???) of residuals
        residStdErr = sqrt(contrastVar); 

        % Convert to T statistics and P values.
        Tstat = contrastFns ./ residStdErr;

        Pval = nan(size(Tstat));
        for j=1:D
            Pval(j,:) = tcdf(Tstat(j,:),dof);
        end
        contrastZ = norminv(Pval);

        t=toc;
        fprintf('Done! Took %.1f seconds.\n',t);


        %% Get group-level stats
        
        [group_RF,group_P] = RunTopLevelGlm_EEG(contrastFns,contrastVar,multcorrect);
        group_Z = norminv(group_P);
        %%
        cd(basedir)
        save(sprintf('%s-%s-%scontrast.mat',experiment,suffix_out,rule),'event_*','contrast*',...
            'group_*','chanlocs','tResponse','experiment','subjects','suffix_in','suffix_out',...
            'rule','sequence','tContrast','legendstr','lambda','multcorrect')



    end

end