function [skeletons, channels] = collectSkeletonData(motionTypes, frac, continuous)

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

skeletons = cell(1,length(subjects));
channels = cell(1,length(subjects));

for i = 1:length(subjects)
    iSub = subjects(i);
    fileNameAsf = sprintf('%02d.asf',iSub);
    if exist(fileNameAsf,'file')
        fprintf('Loading subject %d ...(%d/%d)\n',iSub,i,length(subjects));
        skel = acclaimReadSkel(fileNameAsf);
        skel.subcategory = names{i};
        %         skeletons{i} = acclaimReadSkel(fileNameAsf);
        for j = 1:length(trials{i})
            fileNameAmc = sprintf('%02d_%02d.amc',iSub,j);
            if exist(fileNameAmc,'file')
                numFramesInFile = acclaimNumberOfFrames(fileNameAmc);
                numFramesToExtract = ceil(numFramesInFile*frac);
                if continuous                                                           % Subset of consecutive frame numbers.
                    iSubset = 1:numFramesToExtract;
                else
                    iSubset = sort(randsample(numFramesInFile, numFramesToExtract));    % Subset of uniformly sampled and ordered frame numbers.
                end
                fprintf('\tLoading trial %02d_%02d.amc ...(%d/%d)\n',iSub,j,j,length(trials{i}));
                if frac < 1
                    fprintf('\t\tLoading only %d out of %d frames.\n', numFramesToExtract, numFramesInFile);
                end
                [channel, skeleton] = acclaimSemiLoadChannels(fileNameAmc, skel, iSubset);
                channels{i} = [channels{i} {channel}];
                skeletons{i} = [skeletons{i} {skeleton}];
            end
        end
    end
end

