function [skeletons, channels, xyzChannels] = collectSkeletonData(motionTypes, frac, continuous)

names = {};
subjects = [];
trials = {};

if nargin < 1 || isempty(motionTypes)
    error('Empty or no argument given.')
elseif isstruct(motionTypes)
    motionTypes = {motionTypes};
elseif ~iscell(motionTypes)
    error('Argument must be either a (sub)category struct, or cell array of (sub)category structs. See includeMotionCategories.m for available types.');
end

% Iterate through structs in motionTypes.
for iType = 1:length(motionTypes)
    if ismember('subcategories', fieldnames(motionTypes{iType})) % Category struct.
        for i = 1:length(motionTypes{iType})
            names = [names motionTypes{iType}(i).subcategories.name];
            subjects = [subjects motionTypes{iType}(i).subcategories.subject]; %#ok<*AGROW>
            trials = [trials motionTypes{iType}(i).subcategories.trial];
        end
    elseif ismember('subject', fieldnames(motionTypes{iType})) % Sub-category struct.
        names = [names motionTypes{iType}.name];
        subjects = [subjects motionTypes{iType}.subject];
        trials = [trials motionTypes{iType}.trial];
    end
end

if nargin < 2 || isempty(frac)
    frac = 1;
elseif frac <= 0 || frac > 1
    error('Argument frac must be positive and no larger than 1.')
end

if nargin < 3 || isempty(continuous)
    continuous = true;
end

numSubjects = length(subjects);
skeletons = cell(1, numSubjects);
channels = cell(1, numSubjects);
xyzChannels = cell(1, numSubjects);

for i = 1:numSubjects
    iSub = subjects(i);
    fileNameAsf = sprintf('%02d.asf',iSub);
    if exist(fileNameAsf,'file')
        fprintf('Loading subject %d ...(%d/%d)\n',iSub,i,numSubjects);
        skel = acclaimReadSkel(fileNameAsf);
        skel.subcategory = names{i};
        numTrials = length(trials{i});
        for j = 1:numTrials
            iTrial = trials{i}(j);
            fileNameAmc = sprintf('%02d_%02d.amc',iSub,iTrial);
            if exist(fileNameAmc,'file')
                numFramesInFile = acclaimNumberOfFrames(fileNameAmc);
                numFramesToExtract = ceil(numFramesInFile*frac);
                if continuous                                                           % Subset of consecutive frame numbers.
                    iSubset = 1:numFramesToExtract;
                else
                    iSubset = sort(randsample(numFramesInFile, numFramesToExtract));    % Subset of uniformly sampled and ordered frame numbers.
                end
                fprintf('\tLoading trial %02d_%02d.amc ...(%d/%d)\n',iSub,iTrial,j,numTrials);
                if frac < 1
                    fprintf('\t\tLoading only %d out of %d frames.\n', numFramesToExtract, numFramesInFile);
                end
                [channel, skeleton] = acclaimSemiLoadChannels(fileNameAmc, skel, iSubset);
                skeletons{i} = [skeletons{i} {skeleton}];
                channels{i} = [channels{i} {channel}];
                numFrames = size(channel,1);
                xyzChannel = zeros(numFrames, 93);
                for iFrame = 1:numFrames;
                    xyzChannel(iFrame,:) = reshape( acclaim2xyz(skeleton, channel(iFrame,:)), 1, []);
                end
                xyzChannels{i} = [xyzChannels{i} {xyzChannel}];
            end
        end
    end
end

