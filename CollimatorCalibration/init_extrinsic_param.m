function [rotMatSet, tcp_refine] = init_extrinsic_param(pts3dSet, pts2dSet, nImg, K, tcp)

% Firstly, we estimate the rotation matrix of each image for gobal minimization

for ii = 1:nImg

    pts2d =  pts2dSet{ii};
    pts3d =  pts3dSet{ii};

    % We solve the rotation matrix using the plane homography constraint,
    % or using the dlt-pnp method.
    if 0
        % plane homography constraint
        H = compute_homography(pts2d, pts3d(1:2,:));

        h1 = H(:,1);
        h2 = H(:,2);
        h3 = H(:,3);

        % Eq.(18) in the paper
        lambda1 = 1 / norm(K \ h1);
        lambda2 = 1 / norm(K \ h2);
        lambda = (lambda1 + lambda2) / 2;

        r1 = lambda * ( K \ h1 );
        r2 = lambda * ( K \ h2 );
        r3 = cross(r1, r2);

        R = [r1, r2, r3];

    else
        % dlt-pnp method

        pts3dNorm = pts3d - tcp*ones(1,size(pts3d,2));
        pts2dNorm = inv(K)*[pts2d;ones(1,size(pts2d,2))];
        pts2dNorm = pts2dNorm(1:2,:);

        A = [];
        for jj = 1:size(pts3dNorm,2)
            XYZ = pts3dNorm(:,jj)';
            u = pts2dNorm(1,jj);
            v = pts2dNorm(2,jj);
            A = [A;[XYZ,zeros(1,3),-u.*XYZ;zeros(1,3),XYZ,-v.*XYZ]];
        end

        [~,~,V] = svd(A);
        R = reshape(V(:,end),[3,3])';

    end

    % re-orthogonalize the rotation matrix
    [R_U, R_S, R_V] = svd(R);
    R_init = R_U*R_V';
    if det(R_init)<0
        R_init = -R_init;
    end

    % Refine the rotation matrix to minimize the reprojection error
    param = rodrigues(R_init);

    options = optimset;
    options.Algorithm = 'levenberg-marquardt';
    options.Display = 'none';
    options.MaxIter = 50;
    options.ScaleProblem = 'Jacobian';

    [param_optim] = lsqnonlin(@(param) reproj_error_rot(param, K, tcp, pts3d, pts2d), param, [], [], options);

    R_optim = rodrigues(param_optim);

    rotSet(:,:,ii) = R_optim;
end

% Then, we jointly refine tcp and the rotation matrix of all images by minimizing the reprojection error.

param(1:3,:) = tcp;

for ii = 1:nImg
    rotVec = rodrigues(rotSet(:,:,ii));
    param = [param;rotVec];
end

options = optimset;
options.Algorithm = 'levenberg-marquardt';
options.Display = 'none';
options.MaxIter = 100;
options.ScaleProblem = 'Jacobian';

[param_optim] = lsqnonlin(@(param) reproj_error_rotSet_tcp(param, K, pts3dSet, pts2dSet), param, [], [], options);

tcp_refine = param_optim(1:3);

for ii = 1:nImg
    rotVecOpt = param_optim(3+3*ii-2:3+3*ii);
    rotMatSet(:,:,ii) = rodrigues(rotVecOpt);
end
end

%% Calculate the reprojection error for a single image
function repjError = reproj_error_rot(param, K, tcp, pts3d, pts2d)

R = rodrigues(param);
t = -R*tcp; % obtain from Eq.(4) in the paper
Rt = [R,t;[0,0,0,1]];

pts2dRepj = project_3dto2d(pts3d, Rt, K, [0;0]);
repjError = vecnorm(pts2d - pts2dRepj)'; % reprojection error

end

%% Calculate the reprojection error for all images
function repjErrorSet = reproj_error_rotSet_tcp(param, K, pts3dSet, pts2dSet)

tcp = param(1:3);
nImg = size(pts2dSet,2);

repjErrorSet = [];
for ii = 1:nImg
    pts2d =  pts2dSet{ii};
    pts3d =  pts3dSet{ii};

    rotVec = param(3+3*ii-2:3+3*ii);
    R = rodrigues(rotVec);
    t = -R*tcp;
    Rt = [R,t;[0,0,0,1]];

    pts2dRepj = project_3dto2d(pts3d, Rt, K, [0;0]);

    repjErrorSet = [repjErrorSet;vecnorm(pts2d - pts2dRepj)'];
end
end



