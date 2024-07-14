if size(imgFiles,1) > 0
    for ii = 1:size(pointTxts,1)

        imgPath = strcat(pointTxts(ii).folder, '/' , pointTxts(ii).name(1:end-4), '.bmp');

        if exist(imgPath,"file")
            image = imread(imgPath);
        else
            continue
        end

        if isempty(image)
            continue;
        end

        points2d = points2dSet{ii};

        % show image and 2D points
        figure(ii)
        imshow(image);
        hold on
        plot(points2d(1,:)+1, points2d(2,:)+1, 'r.', 'MarkerSize', 10);
        title(['I' pointTxts(ii).name(2:end-4) ' and 2D feature points']);
        axis on
        hold off
        drawnow
    end
end