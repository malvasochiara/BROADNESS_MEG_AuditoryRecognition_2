function[OUT] = clusterbased_permutationtest(S)
% Performs a T-test on original data and identifies the clusters. Then
% performs a cluster- based permutation test and saves the results.
%%%%%%%%%%%%%%%%%
%  INPUT : - S.reference = time series of the reference. !!!! Dimensions are time x  subjects !!!!
%          - S.testing   = time series you want to test against reference. !!!! Dimensions are time x  subjects !!!!
%          - S.nperm     = number of permutations
%          - S.alpha     = significance level
%          - S.threshold = significance level you consider significant after permutation test
%          - S.time      = array containing the time points
%          - S.stattype  = how to compute the overall statistics of each cluster: .'max' : (default) maximum (absolute) tvalue
%                                                                                 .'size': size of the cluster
%                                                                                 .'sum' : sum of the tvalues
%
% OUTPUT : -OUT.dataastat = T values of your initial data
%          -OUT.nullstat  = T values of all the permutations
%


% Developed by Chiara Malvaso, chiara.malvaso@studio.unibo.it
% Supervised by Leonardo Bonetti, leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk 


    reference = S.reference;
    testing = S.testing;
    nperm = S.nperm;
    alpha = S.alpha;
    threshold = S.threshold;
    time = S.time;
    if isempty(S.stattype)
        stattype = 'max';
    else
    stattype = S.stattype;
    end
    %PCs = size(dataa,1);
    timepoints = size(reference,1);
    %conditions = size(dataa,3);
    %comparisons = conditions -1; % number of comparisons you make
    subjects = size(reference,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1: STATISTIC ON ORIGINAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the T test for each time point
    T_storing = zeros(timepoints,1);
    P_storing = zeros(timepoints,1);
    for tt = 1:timepoints
        [~,p,~,stats] = ttest(squeeze(reference(tt,:)),squeeze(testing(tt,:)));
        P_storing(tt) = p;
        T_storing(tt) = stats.tstat;
    end
% Identify clusters on binary data
    cluster_originaldata  = []; %contains informations about the index, size etc of clusters
    cluster_originaldata  = bwconncomp(double(P_storing < alpha));
    
% Storing the size and the maximum T value of each cluster WITH THE SIGN
    cluster_size = zeros(cluster_originaldata.NumObjects,1);
    cluster_max = zeros(cluster_originaldata.NumObjects,1);
    cluster_stats_originaldata = zeros(cluster_originaldata.NumObjects,1); %store the statistic of each cluster 
    
    for cc = 1:cluster_originaldata.NumObjects % find the overall statistic for each cluster so you can compare it to the distribution you'll get from section 2
        index = cluster_originaldata.PixelIdxList{cc};
        cluster_size(cc) = length(index); %storing the size
        [absmax, indexmax] = max(abs(T_storing(index)));
        cluster_max(cc) = T_storing(index(1)-1+indexmax); %storing the Tvalue with the sign
        if strcmp(stattype,'max')
            cluster_stats_originaldata(cc) = absmax;
        else
            if strcmp(stattype,'size')
                cluster_stats_originaldata(cc) = cluster_size(cc);
            else
                cluster_stats_originaldata(cc) = sum(abs(T_storing(index)));
            end
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2: PERMUTATION TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T_cluster = zeros(nperm,1); % store one T value for each permutation
    for n = 1:nperm % iterate the testing and cluster identification process for each permutation
        %disp(n)
        T_storing = zeros(timepoints,1);
        P_storing = zeros(timepoints,1);
        for tt = 1:timepoints  %repeate the permutation and the T test for each time point
            perm_reference = zeros(subjects,1);
            perm_testing = zeros(subjects,1);

            for sub = 1:subjects %compute the permutation for each subject
                if rand > 0.5
                    perm_reference(sub) = reference(tt, sub);
                    perm_testing(sub) = testing(tt, sub);
                else
                    perm_reference(sub) = testing(tt, sub);
                    perm_testing(sub)= reference(tt, sub);
                end

            end
            % perform the t test for each time point and store the output
            [~,p,~,stats] = ttest(squeeze(perm_reference),squeeze(perm_testing));
            P_storing(tt) = p;
            T_storing(tt) = stats.tstat;

        end
        clear cluster
        cluster = bwconncomp(double(P_storing < alpha));
        
        clear cluster_stats
        for cc = 1:cluster.NumObjects % find the overall statistic for each cluster
            index = cluster.PixelIdxList{cc};
            if strcmp(stattype,'max')
                cluster_stats(cc) = max(abs(T_storing(index)));
            else
                if strcmp(stattype,'size')
                    cluster_stats(cc) = length(index);
                else
                    cluster_stats(cc) = sum(abs(T_storing(index)));
                end
            end
        end
        % Storing the maximum statististic recorder for each cluster
        T_cluster(n,1) = max(cluster_stats);


    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 3: COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the adjusted p values for each cluster in the original data:
    OUT = {};
    sig_number = 0;
    for cc = 1:length(cluster_stats_originaldata) %going over the clusters identified in the original data
        count = 0;
        for n = 1:nperm
            if T_cluster(n) >= cluster_stats_originaldata(cc)
                count = count + 1;
            end
        end
        p_val = count/nperm;
        % check if the pvalue is significant
        
        if p_val < threshold
            sig_number = sig_number +1;
            OUT(sig_number).clustersize = cluster_size(cc);
            OUT(sig_number).pvalue = p_val;
            timeindex = cluster_originaldata.PixelIdxList{cc};
            OUT(sig_number).temporalextent = time(timeindex(end)) - time(timeindex(1));
            OUT(sig_number).time = [time(timeindex(1)) time(timeindex(end))];
%           OUT(sig_number).stats = cluster_stats_originaldata(cc);
            OUT(sig_number).Tvalue = cluster_max(cc);

            
        end
    end
    
    



end
