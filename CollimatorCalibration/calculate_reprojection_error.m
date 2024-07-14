function [meanError, medianError,stdError, repjErrorSet] = calculate_reprojection_error(K, poseSet, d, pts3dSet, pts2dSet, nImg)

repjErrorSet = [];

if nImg == 1
    pts2dSet{1} = pts2dSet;
    pts3dSet{1} = pts3dSet;
end


for ii = 1:nImg

    pts2d =  pts2dSet{ii};
    pts3d =  pts3dSet{ii};

    Rt = poseSet(:,:,ii);

    pts2dRepj = project_3dto2d(pts3d, Rt, K, d);

    repjError = vecnorm(pts2d - pts2dRepj);

    repjErrorSet = [repjErrorSet;repjError'];
end

meanError = mean(repjErrorSet);
medianError = median(repjErrorSet);
stdError = std(repjErrorSet);

end