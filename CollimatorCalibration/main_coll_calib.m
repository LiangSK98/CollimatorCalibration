function calibResult = main_coll_calib(points3dSet, points2dSet, config)

nImg = config.nImg;

if(size(points3dSet,2) ~= nImg || size(points2dSet,2) ~= nImg)
    repjError('The dimension of 3D/2D point set is not equal to the number of images.\n');
end

if (isfield(config,'useCenter') && config.useCenter == 1)
    cxcy = config.imageSize ./ 2;
else
    cxcy = [];
end

% Initialization of the intrinsic parameters:
if nImg > 2
    [K_init, tcp_init] = init_intrinsic_param(points3dSet, points2dSet, nImg, cxcy);
elseif nImg == 2
    [K_init, tcp_init] = init_intrinsic_param_two_view(points3dSet, points2dSet, nImg, cxcy);
else
    repjError('calibration requires at least 2 iamges!')
end

d_init = [0;0];
fprintf('- Initialization two order distortion: [ %3.4f  %3.4f ]\n', d_init)

% Initialization of poses for all calibration images
[rotMatSet_init, tcp_init] = init_extrinsic_param(points3dSet, points2dSet, nImg, K_init, tcp_init);

% report the initial estimation of the camera parameters

fprintf('\nInitial guess of camera parameters: \n');
fprintf('- Focal Length:                [fx, fy] = [ %3.4f, %3.4f ]\n', K_init(1,1), K_init(2,2))
fprintf('- Principal point:             [cx, cy] = [ %3.4f, %3.4f ]\n', K_init(1,3), K_init(2,3))
fprintf('- Distortion:                  [d1, d2] = [ %3.4f, %3.4f ]\n', d_init)
fprintf('- Fixed camera position:            tcp = [ %3.4f, %3.4f, %3.4f ]\n', tcp_init)

% ----------- Main refinement steps --------------

% Initialization of the parameter vector

fprintf('\nWaiting for the calibration refinement step ......\n')

param_init(1,:) = K_init(1,1); % fx
param_init(2,:) = K_init(2,2); % fy
param_init(3,:) = K_init(1,3); % cx
param_init(4,:) = K_init(2,3); % cy
param_init(5:6,:) = d_init;    % distoration
param_init(7:9,:) = tcp_init;  % tcp

for ii = 1:nImg
    rotMat = rotMatSet_init(:,:,ii);
    rotVec = rodrigues(rotMat);
    param_init = [param_init;rotVec];
end

options = optimset;
options.Algorithm = 'levenberg-marquardt';
options.Display = 'none';
options.TolFun = 1e-8;
options.TolX = 1e-8;
options.MaxIter = 100;
options.ScaleProblem = 'Jacobian';

[param_optim] = lsqnonlin(@(param) reprojection_error_BA(param, points3dSet, points2dSet, nImg), param_init, [], [], options);

% Extract the refined camera parameters

fx = param_optim(1);
fy = param_optim(2);
cx = param_optim(3);
cy = param_optim(4);

K_optim = [fx, 0, cx;
    0,  fy,cy;
    0,  0, 1];

d_optim = [param_optim(5);param_optim(6)];
tcp_optim = param_optim(7:9);

for ii = 1:nImg
    rotVec = param_optim(9+3*ii-2:9+3*ii);
    rotMat = rodrigues(rotVec);
    rotMatSet_optim(:,:,ii) = rotMat;
    R = rotMat;
    t = -R*tcp_optim;
    poseSet(:,:,ii) = [R,t;[0,0,0,1]];
end

% Calculate reprojection error
[~,~,~,repjError] = calculate_reprojection_error(K_optim, poseSet, d_optim, points3dSet, points2dSet, nImg);

calibResult.K = K_optim;
calibResult.d = d_optim;
calibResult.tcp = tcp_optim;
calibResult.rotMatSet = rotMatSet_optim;
calibResult.poseSet = poseSet;
calibResult.repjError = repjError;

% When a low-precision collimator is used, refining the position for each
% image is useful to precisely update the camera parameters.

if (isfield(config,'optimPos') && config.optimPos == 1)

    param_init2(1,:) = K_optim(1,1); % fx
    param_init2(2,:) = K_optim(2,2); % fy
    param_init2(3,:) = K_optim(1,3); % cx
    param_init2(4,:) = K_optim(2,3); % cy
    param_init2(5:6,:) = d_optim;    % distoration

    for ii = 1:nImg
        rotMat = rotMatSet_optim(:,:,ii);
        rotVec = rodrigues(rotMat);
        param_init2 = [param_init2;rotVec;tcp_optim];
    end

    options = optimset;
    options.Algorithm = 'levenberg-marquardt';
    options.Display = 'none';
    options.TolFun = 1e-8;
    options.TolX = 1e-8;
    options.MaxIter = 100;
    options.ScaleProblem = 'Jacobian';

    [param_optimPos] = lsqnonlin(@(param) reprojection_error_optimPos(param, points3dSet, points2dSet, nImg), param_init2, [], [], options);

    fx = param_optimPos(1);
    fy = param_optimPos(2);
    cx = param_optimPos(3);
    cy = param_optimPos(4);

    K_optimPos = [fx, 0, cx;
        0,  fy,cy;
        0,  0, 1];

    d_optimPos = [param_optimPos(5);param_optimPos(6)];

    for ii = 1:nImg
        rotVec = param_optimPos(6+6*ii-5:6+6*ii-3);
        rotMat = rodrigues(rotVec);
        rotMatSet_optimPos(:,:,ii) = rotMat;
        tcp = param_optimPos(6+6*ii-2:6+6*ii);

        R = rodrigues(rotVec);
        t = -R*tcp;
        poseSet(:,:,ii) = [R,t;[0,0,0,1]];
    end

    % Calculate reprojection error
    [~,~,~,repjError] = calculate_reprojection_error(K_optimPos, poseSet, d_optimPos, points3dSet, points2dSet, nImg);

    calibResult.K = K_optimPos;
    calibResult.d = d_optimPos;
    calibResult.tcp = tcp_optim;
    calibResult.rotMatSet = rotMatSet_optimPos;
    calibResult.poseSet = poseSet;
    calibResult.repjError = repjError;
end

fprintf('\nCamera parameters after refinement steps: \n')
fprintf('- Focal Length:                [fx, fy] = [ %3.4f, %3.4f ]\n', calibResult.K(1,1), calibResult.K(2,2))
fprintf('- Principal point:             [cx, cy] = [ %3.4f, %3.4f ]\n', calibResult.K(1,3), calibResult.K(2,3))
fprintf('- Distortion:                  [d1, d2] = [ %3.4f, %3.4f ]\n', calibResult.d)
fprintf('- Fixed camera position:            tcp = [ %3.4f, %3.4f, %3.4f ]\n', calibResult.tcp)
fprintf('- Reprojection error:  [mean, med, std] = [ %3.5f, %3.5f, %3.5f ]\n', ...
    mean(calibResult.repjError), median(calibResult.repjError), std(calibResult.repjError));
end