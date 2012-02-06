% DEMCMUMOCAPEMRCA1 EM-RCA demo on reconstruction of the stick man from 3D
% sensor data from motions across the CMU mocap database.
%
% FORMAT
% DESC
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2011, 2012
%
% RCA

clc, clear
addpath(genpath('~/mlprojects/matlab/general/'))
addpath(genpath('~/mlprojects/rca/matlab/glasso/'))
addpath(genpath('~/mlprojects/rca/matlab/L1General'))
importTool({'rca','ndlutil','datasets','mocap'})
addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
includeMotionCategories

[walkSkeletons, walkChannels] = collectSkeletonData(walking, .1, true);
[runSkeletons, runChannels] = collectSkeletonData(running, .1, true);
[jumpSkeletons, jumpChannels] = collectSkeletonData(jumping, .1, true);
[miscSkeletons, miscChannels] = collectSkeletonData({playground, physical_activities_and_sports, situations_and_scenarios }, .1, true);

% [danceChannels, danceSkeletons] = acclaimSemiLoadChannels('02_05.amc', acclaimReadSkel('02.asf'));
[danceSkeletons, danceChannels] = collectSkeletonData(dance, .05, false);
skelPlayData(danceSkeletons{1}{1}, danceChannels{1}{1});

% Animate all extracted motions.
for i = 1:length(danceSkeletons)
    for j = 1:length(danceSkeletons{i})
        skelPlayData(danceSkeletons{i}{j}, danceChannels{i}{j});
    end
end

