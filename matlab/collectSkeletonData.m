% function collectSkeletonData()

addpath(genpath('~/Desktop/CMUmocap/all_asfamc/subjects/'))
connect = cell(1,143);

sumConnect = zeros(31,31);
for i = 1:143
    if i <= 9
        fileNameAsf = ['0' num2str(i) '.asf'];
    else
        fileNameAsf = [num2str(i) '.asf'];
    end
    
    if exist(fileNameAsf,'file')
        skel = acclaimReadSkel(fileNameAsf);
        connect{i} = skelConnectionMatrix(skel);
        sumConnect = sumConnect + connect{i};
        imagesc(sumConnect), colorbar;
        [channels, skel] = acclaimLoadChannels(fileNameAmc, skel);
    end
end