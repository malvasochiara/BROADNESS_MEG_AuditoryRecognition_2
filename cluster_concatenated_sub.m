function OUT = cluster_concatenated_sub(S)
% Load the necessary data and concatenate subjects, then perform pca on the 
%  concatenated matrix and reconstruct the time serie for each subject.
% Variance and time series are saved in a selected folder.
% INPUT:
%         -S structure with fields:
%           -outputdir  = path to the folder were variance and time series are going to be saved
%           -folderpath = path to the folder containing the data for each subject
%           -sign_eig   = method for normalizing eigenvectors sign (SUGGESTED EITHER 'max_abs' or 'average'):
%                            -'occurrences' = using mean of the negative/positive values occurrences
%                            -'max_abs' = on the basis of the sign of the maximum value in absolute terms
%                            -'average' = on the basis of the sign of the average of weights
% OUTPUT:
%        - OUT structure with fields:
%           -W             = weight coeffiecients (normalized sign)
%           -variance_PCS  = explained variance


% Developed by Chiara Malvaso, chiara.malvaso@studio.unibo.it
% Supervised by Leonardo Bonetti, leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk 

OUT = [];

% extracting inputs
outputdir = S.outputdir;
folderpath = S.folderpath;


pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl); % add to path LBPD and other OSL functions
addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/Broadness')
%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

%loading data in 3559-voxel space (8mm)
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/sources_main_effects.mat');

%loading single subjects
fileList = dir([folderpath '/SUBJ*.mat']);
N = length(fileList);
bigMatrix = zeros (3559,1026,N); % only consider condition 1
% 
% for i = 1:length(fileList)
for i = 1:N
    disp(['Loading subject ' num2str(i)]);
    currentFileName = fileList(i).name;
    currentFilePath = fullfile(folderpath, currentFileName);
    try
        data = load(currentFilePath);
        bigMatrix(:,:,i) = squeeze(data.OUT.sources_ERFs(:,:,1)); %selcting condition 1
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            disp('errore in %s', currentFilename);
            continue
        end
    end
end
%actual computation
%adjusting polarity
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end

dum = zeros(size(bigMatrix,1),size(bigMatrix,2),N);
for jj = 1:size(t_val_s,1) %over brain voxels
    for s = 1:N % over subjects
        dum(jj,:,s) = bigMatrix(jj,:,s) .* vect(jj,1); %reversing (or not)..
        disp([' - source ' num2str(jj)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H =  reshape(dum, size(dum,1), size(dum,2)*size(dum,3));
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
% normalizing wcoeff by multiplying for the covariance matrix
C = cov(H'); 

OUT.W_normalized = (VV2')*C; %weights of PCA
OUT.variance_PCS = vare; %storing variance of significant PCs

save([outputdir '/pca_concatenatedsubjects_COND1_varexp'], 'vare');
disp ('Variance was saved');

PCs = 1:2; %number of components you want to see


J = zeros(1025,PCs(end),N);

    for ii = 1:N %over subjects
        J(:,:,ii) = dum(:,1:1025,ii)' * VV2(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
        disp(['Computing time serie ' num2str(ii)]);
    end
    
% saving:
save([outputdir '/pca_concatenatedsubjects_timeserie'], 'J');
disp('Time series have been saved')
end



