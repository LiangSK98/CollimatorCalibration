function repjErrorSet = reprojection_error_BA(param, pts3dSet, pts2dSet, nImg)

% Unpacking parameter vector

fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);

K = [fx, 0, cx;
     0,  fy,cy;
     0,  0, 1];

d = [param(5);param(6)];
tcp = param(7:9);

% Calculate cost function
repjErrorSet = [];
for ii = 1:nImg

    pts2d =  pts2dSet{ii};
    pts3d =  pts3dSet{ii};

    rotVec = param(9+3*ii-2:9+3*ii);
    R = rodrigues(rotVec);
    t = -R*tcp;
    Rt = [R,t;[0,0,0,1]];

    pts2dRepj = project_3dto2d(pts3d, Rt, K, d);
    
    repjError = vecnorm(pts2d - pts2dRepj);
    repjErrorSet = [repjErrorSet;repjError'];
end
% Robust Cauchy cost function

c=1.385;
repjErrorSet=sqrt((c^2/2)*log(1+(repjErrorSet/c).^2));
end