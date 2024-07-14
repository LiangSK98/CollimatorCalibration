
% Tested demo code on MATLAB R2021b

clear all; clc; close all;

addpath CollimatorCalibration

fprintf('\nHere is an example of camera calibration using collimator system! \n\n')

% Load images (.bmp) and text files containing 3D and 2D points (.txt )

dataPath = './data/data2/'; % or './data/data2/'
imgFiles = dir([dataPath '*.bmp']);
pointTxts = dir([dataPath '*.txt']);

% Random number and random order (Users can also be commented out.)

% index = randperm(20);
% pointTxts = pointTxts(index(1:randi([5,20])),:);

for ii = 1:size(pointTxts,1)

    % load 3D and 2D points
    points2d3d = load(strcat(dataPath,pointTxts(ii).name));

    points2d = points2d3d(:,1:2)';
    points3d = points2d3d(:,3:4)';
    points3d(3,:) = zeros(1,size(points3d,2));

    points2dSet{ii} = points2d;
    points3dSet{ii} = points3d;

    % load images (if have)

    if ~isempty(imgFiles)
        imgPath = strcat(dataPath, pointTxts(ii).name(1:end-4), '.bmp');
        if exist(imgPath,'file')
            image = imread(imgPath);
            imageSize(ii,:) = [size(image,2), size(image,1)];
        else
            imageSize(ii,:) = nan;
        end
    end
end

% Show image and 2D feature points (Not necessary, can be commented out)

if ~isempty(imgFiles)

%     show_data()
%     fprintf('\nPress any key to continue ...... \n\n')
%     pause

    close all
end

% Calibration configuration

% Check image size
if exist('imageSize')
    if ((max(imageSize(:,1)) - min(imageSize(:,1)))~=0 & (max(imageSize(:,2)) - min(imageSize(:,2)))~=0)
        error('All images should have the same size!')
    end
    config.imageSize = imageSize(1,:); % for /data1
else
    config.imageSize = [1080,960]; % for /data2
end

% Number of images for calibration
config.nImg = size(pointTxts,1);

% Whether to use the image center as the initial principal point
config.useCenter = 1;

% Whether to refine the position of each image
config.optimPos = 1;

fprintf('Calibration configuration: \n')
fprintf('- Number of images: %d\n' ,config.nImg)
fprintf('- Resolution of the image: [ %d  %d ] \n', config.imageSize)
fprintf('- Initialization of the principal point at the center of the image?  ')
if (isfield(config,'useCenter') && config.useCenter == 1)
    fprintf('yes\n')
else
    fprintf('no\n')
end
fprintf('- Refine the position for each image?  ')
if(isfield(config,'optimPos') && config.optimPos == 1)
    fprintf('yes\n')
else
    fprintf('no\n')
end

% Run the main calibration routine
calibResult = main_coll_calib(points3dSet, points2dSet, config);

% Shows the distribution of reprojection errors
show_reprojection_error()

% Thanks
