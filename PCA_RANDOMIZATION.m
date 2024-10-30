function [ OUT ] = PCA_RANDOMIZATION(S)
% Load the data for each subject and perform the selected number of
% permutations. Pca is performed on randomized data and the corresponding
% time series are computed. The explained variance and the time series on
% randomized data are saved in the selected folder. Pca is also computed on
% main effect and the explained variance is saved.
%
%
% INPUT:
%       -S structure with fields:
%           -maineffect   = matrix that is going to be used to compute pca (brain sources, time points, condtions)
%           -folderpath   = path to the folder containing the data for each subject
%           -rand_l       = strategy for randomisation:
%                           1 = randomising only time-points (independently for each time series)
%                           2 = randomising only space (independently for each time-point) 
%                           3 = randomising both time and space
%                           4 = randomising the spatial labels while keeping the order of time points
%           -permnum      = number of permutations
%           -outputdir    = path to the folder where outputs must be saved
%           -vect         = vector to adjust the polarity computed on main effect
%
% OUTPUT:
%        - OUT structure with fields:
%           -wcoeff                   = weight coeffiecients from pca on main effect
%           -timeserie_rand_averaged  = time series recontructed on randomized data
%           -wcoeff_rand              = weight coefficient from pca on randomized data


% Developed by Chiara Malvaso, chiara.malvaso@studio.studio.unibo.it
% Supervised by Leonardo Bonetti, leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk 
%%

%initializing output structure
    OUT = [];
    %extracting inputs
    outputdir = S.outputdir;
    folderpath = S.folderpath;
    permnum = S.permnum;
    if (permnum <= 0)
        error('Select a positive number of permutations')
    end
    rand_l = S.rand_l;
    if ~ismember(rand_l, [1 2 3 4])
        error('Select an existent permutation strategy')
    end
    H = mean(S.maineffect(:,:,:),3);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA COMPUTATION ON ORIGINAL DATA + CORRECTION OF WCOEFFS
    [wcoeff,~,~,~,vare] = pca(H');
    %normalizing eigenvectors signs
    dumones = ones(size(wcoeff,1),size(wcoeff,2)); %vector of 1s with lenngth of significant PCs
    switch S.sign_eig
        case 'occurrences'
            % 1) normalizing eigenvectors sign by using mean of the negative/positive values occurrences
            dumones(:,mean(wcoeff > 0) < 0.5) = -1; %assigning -1 to eigenvectors that have more positive than negative weights
        case 'max_abs'
            % 2) normalizing eigenvectors sign on the basis of the sign of the maximum value in absolute terms
            [~,mi] = max(abs(wcoeff)); %getting maximum values indices
            ab = zeros(1,length(mi));
            for ii = 1:length(mi) %over significant PCs
                ab(1,ii) = wcoeff(mi(ii),ii); %storing original values with signs (corresponding to maximum in absolute terms)
            end
            dumones(:,sign(ab)<0) = -1;
        case 'average'
            % 3) normalizing eigenvectors sign on the basis of the sign of the average of weights
            mVV = mean(wcoeff);
            dumones(:,sign(mVV)<0) = -1;
    end
    VV2 = wcoeff .* dumones; %getting proper sign for eigenvectors
    
    % storing data
    OUT.wcoeff = VV2;

%%%%%%%%%%%%%%%%%%%%%%% loading single subject data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loading data in 3559-voxel space (8mm)
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/sources_main_effects.mat');

fileList = dir([folderpath '/SUBJ*.mat']);
N =length(fileList);
bigMatrix = zeros (3559,1026,N, 5);

for i = 1:N
    disp(['Loading subject ' num2str(i)]);
    currentFileName = fileList(i).name;
    currentFilePath = fullfile(folderpath, currentFileName);
    try
        for cc = 1:5
        data = load(currentFilePath);
        bigMatrix(:,:,i, cc) = squeeze(data.OUT.sources_ERFs(:,:,cc));
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            disp('errore in %s', currentFilename);
            continue
        end
    end
end

%actual computation
%adjusting polarity
vect = S.vect;

H_sub = zeros(size(bigMatrix,1),size(bigMatrix,2),N,5);
for jj = 1:size(t_val_s,1) %over brain voxels
    for s = 1:N % over subjects
        for cc = 1:5 %over conditions
            H_sub(jj,:,s,cc) = bigMatrix(jj,:,s,cc) .* vect(jj,1); %reversing (or not)..
            disp([' - source ' num2str(jj)])
        end
    end
end   
    
    
    
%%%%%%%%%%%%%%%%%%%%%%% PCA on random data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_rand = zeros(size(wcoeff,2),2,size(H_sub,3), size(H_sub,4),permnum); % 1025 x 2 x subjects x conditions x permutations
    varef = zeros(size(H,2)-1,permnum); % 1025 x permutations
    wcoeff_rand = zeros(size(H,1), size(H,2)-1,permnum); % 3559 x 1025 x permutations
    for pp = 1:permnum %over permutations
        disp(['Permutation ' num2str(pp)])
        clear r_resh;
        clear r_resh_ts;
        if rand_l == 1     %randomizing only time
            r_resh_ts = zeros(size(H,1),size(H,2),size(H_sub,3), size(H_sub,4));
            r_resh = zeros(size(H,1), size(H,2));
            for ss = 1:size(H,1) %over brain sources
                %rng('shuffle');
                clear idx_dummy;
                idx_dummy = randperm(size(H,2)); %create a permuted array of the indices of the temporal dimension
                r_resh(ss,:) = H(ss, idx_dummy);
                for sub = 1:size(H_sub,3) %over subjects
                    for cc = 1:size(H_sub,4) %over conditions
                        r_resh_ts(ss,:,sub,cc) = H_sub(ss,idx_dummy,sub,cc);
                    end
                end
            end
        elseif rand_l == 2 %randomizing only space
                r_resh_ts = zeros(size(H,1),size(H,2),size(H_sub,3), size(H_sub,4));
                r_resh = zeros(size(H,1), size(H,2));
             for ss = 1:size(H,2) %over time-points
                %rng('shuffle');
                clear idx_dummy;
                idx_dummy = randperm(size(H,1)); %create a permuted array of the indices of the spatial dimension
                r_resh(:,ss) = H(idx_dummy,ss);
                    for sub = 1:size(H_sub,3) %over subjects
                        for cc = 1:size(H_sub,4) %over conditions
                            r_resh_ts(:,ss,sub,cc) = H_sub(idx_dummy,ss,sub,cc);
                        end
                    end
            end
        elseif rand_l == 3 %randomizing both space and time
            %randomizing data..
            r_resh_ts = zeros(size(H_ts_rand,1),size(H_ts_rand,2),size(H_ts_rand,3));
            idx_dummy = randperm(size(H,1)*size(H,2)); %create a permuted array from 1 to size of original data vector
            r_dummy = zeros(size(H,1)*size(H,2),1); %initialise new vector
            r_dummy(1:size(H,1)*size(H,2)) = H(idx_dummy); %taking element in matrix data with index idx_dummy

            for ii = 1:size(H_ts_rand,3)
                r_resh_ts(:,:,ii) = reshape(r_dummy,[size(H,1),size(H,2)]); %reshaping the vector into a matrix shaped as the original data
            end
            %running PCA in single steps on randomised data
        elseif rand_l == 4
                r_resh_ts = zeros(size(H,1),size(H,2),size(H_sub,3), size(H_sub,4));
                r_resh = zeros(size(H,1), size(H,2));
                index = randperm(size(H,1));
            for jj = 1:length(index)
                r_resh(jj,:) = H(index(jj),:);
                for sub = 1:size(H_sub,3) %over subjects
                        for cc = 1:size(H_sub,4) %over conditions
                            r_resh_ts(jj,:,sub,cc) = H_sub(index(jj),:,sub,cc);
                        end
                end
            end
        end
        %r_resh_averaged = mean(r_resh(:,:,:),3);
        clear wcoeff_temp;
        clear varef_temp;
        [wcoeff_temp,~,~,~,varef_temp] = pca(r_resh');
        varef(:,pp) = varef_temp;
        wcoeff_rand(:,:,pp) = wcoeff_temp;
        disp(['perm ' num2str(pp) ' varef ' num2str(max(varef_temp))]);
        for sub = 1:size(H_sub,3) % over subjects
            for cc = 1:size(H_sub,4) % over conditions 
                J_rand(:,:,sub, cc,pp) =(r_resh_ts(:,1:size(wcoeff_rand,2),sub, cc)' * wcoeff_temp(:,1:2)); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
            end
        end        
    end
    varef_averaged = mean(varef(:,:),2); %averaging over permutations
    wcoeff_rand = mean(wcoeff_rand(:,:,:),3); %averaging over permutations
    J_rand_averaged = mean(J_rand(:,:,:,:,:),5); %averaging over permutations

    %storing outputs from random data
    OUT.timeserie_rand_averaged = J_rand_averaged;
    OUT.W_RAND = wcoeff_rand;

%%%%%%%%%%%%% saving %%%%%%%%%%%%%%%%%%%%%%%%%
if rand_l == 1
    save([outputdir '/varexp_from_maineffect'], 'vare');
    save([outputdir '/randomiazionoverTIME_varexp'], 'varef_averaged');
    save([outputdir '/randomiazionoverTIME_tsrand_eachperm'], 'J_rand');
else
    if rand_l == 2
        save([outputdir '/varexp_from_maineffect'], 'vare');
        save([outputdir '/randomiazionoverSPACE_varexp'], 'varef_averaged');
        save([outputdir '/randomiazionoverSPACE_tsrand_eachperm'], 'J_rand');
    else
        if rand_l == 3
            save([outputdir '/varexp_from_maineffect'], 'vare');
            save([outputdir '/randomiazionoverSPACEANDTIME_varexp'], 'varef_averaged');
            save([outputdir '/randomiazionoverSPACEANDTIME_tsrand_eachperm'], 'J_rand');
            
        else
            
            save([outputdir '/varexp_from_maineffect'], 'vare');
            save([outputdir '/randomiazionoverSPACELABELS_varexp'], 'varef_averaged');
            save([outputdir '/randomiazionoverSPACELABELS_tsrand_eachperm'], 'J_rand');
        end
    end

end