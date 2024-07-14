close all;
% show the reprojection error of the calibraton result

if (size(imgFiles,1) > 0)
    noImg = 0;
else
    noImg = 1;
end

% load camera parameters
K = calibResult.K;
d = calibResult.d;

for ii = 1:config.nImg

    % Recalculate reprojection error
    pts2d = points2dSet{ii};
    pts3d = points3dSet{ii};

    Rt = calibResult.poseSet(:,:,ii);

    pts2dRepj = project_3dto2d(pts3d, Rt, K, d);

    repjErrorSet{ii} = pts2dRepj - pts2d;
    meanRepjError(ii) = mean(vecnorm(pts2dRepj - pts2d));

    % show image (if have) and image points
    figure(ii)
    if ~noImg
        imgPath = strcat(pointTxts(ii).folder, '/' , pointTxts(ii).name(1:end-4), '.bmp');
        image = imread(imgPath);
        imshow(image);
        hold on
    end
    plot(pts2d(1,:), pts2d(2,:),'go');
    hold on
    plot(pts2dRepj(1,:), pts2dRepj(2,:), 'r+');
    title(['I' pointTxts(ii).name(2:end-4) '- image points (o) and reprojected points (+)'])
    axis on
    axis([1 config.imageSize(1) 1 config.imageSize(2)])
    hold off;
    drawnow
end

% Show the distribution of reprojection errors
colors = 'brgkcm';

figure(config.nImg+1)
for ii = 1:config.nImg
    repjError = repjErrorSet{ii};
    plot(repjError(1,:), repjError(2,:), [colors(rem(ii-1,6)+1), '+'])
    hold on
end
hold off
axis equal
title('Distribution of reprojection error ( pixel )')
xlabel('x');
ylabel('y');
drawnow;

figure(config.nImg+2)
bar(meanRepjError)
title('Mean reprojection error of each image')
xlabel('Images')
ylabel('Mean reprojection error (pixel)')
set(gca,'XTick',[1:config.nImg])