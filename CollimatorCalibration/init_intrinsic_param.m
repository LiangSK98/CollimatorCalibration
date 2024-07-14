function [K, tcp] = init_intrinsic_param(points3dSet, points2dSet, nImg, cxcy)

% Initialize the homographies

homoSet = zeros(3,3,nImg);
invHomoSet = zeros(3,3,nImg); % The inverse of homography matrix
lambda_i1 = ones(1,nImg); % lambda_i1 can be computed by Eq.(10) in the paper

for ii = 1:nImg

    pts2d =  points2dSet{ii};
    pts3d =  points3dSet{ii};

    H = compute_homography(pts2d, pts3d(1:2,:));
    homoSet(:,:,ii) = H;
    invHomoSet(:,:,ii) = inv(H);

    % Calculate lambda_i1
    H1 = homoSet(:,:,1);
    lambda_i1(ii) = det(H1\H)^(1/3);
end

% Do not initialization of the principal point at the center of the image.
if isempty(cxcy)

    % Coefficient matrics

    V = []; % V = [v11,v12,v13,v22,v23,v33];
    b = []; % b = [b11,b12,b13,b22,b23,b33];

    % Eq.(15) in the paper
    vij = @(i,j,H) [ H(i,1)*H(j,1);
        H(i,1)*H(j,2) + H(i,2)*H(j,1);
        H(i,1)*H(j,3) + H(i,3)*H(j,1);
        H(i,2)*H(j,2);
        H(i,2)*H(j,3) + H(i,3)*H(j,2)]';
    bij = @(i,j,H)[- H(i,3)*H(j,3)];

    for ii = 1:nImg

        invH = invHomoSet(:,:,ii);

        v11 = vij(1, 1, invH);
        v12 = vij(1, 2, invH);
        v13 = vij(1, 3, invH);
        v22 = vij(2, 2, invH);
        v23 = vij(2, 3, invH);
        v33 = vij(3, 3, invH);

        b11 = bij(1, 1, invH);
        b12 = bij(1, 2, invH);
        b13 = bij(1, 3, invH);
        b22 = bij(2, 2, invH);
        b23 = bij(2, 3, invH);
        b33 = bij(3, 3, invH);

        V_ii = [[v11;v12;v13;v22;v23;v33],-1/(lambda_i1(ii)^2) * eye(6)];
        b_ii = [b11;b12;b13;b22;b23;b33];

        V = [V;V_ii];
        b = [b;b_ii];
    end

    % solving linear equations
    wa = V\b; % The least square solution of V*wa = b. Same as inv(V'*V)*V'*b
    w = wa(1:5);
    a = wa(6:11);

    %  Image of the Absolute Conic
    W = [w(1), w(2), w(3);
        w(2), w(4), w(5);
        w(3), w(5), 1];

    cx = W(1,3);
    cy = W(2,3);
    fy = sqrt(W(2,2) - cy*cy);
    fx = sqrt(W(1,1) - cx*cx);

    % initial guess of intrinsic

    K = [fx, 0, cx;
        0, fy, cy;
        0, 0, 1];

    % The position of the camera in the calibration target coordinate system
    A = [a(1), a(2), a(3);
        a(2), a(4), a(5);
        a(3), a(5), a(6)];
    A = A ./ A(3,3);

    % tcp
    tcp_x = A(1,3);
    tcp_y = A(2,3);
    tcp_r = sqrt(A(1,1) - A(1,3)^2);
    tcp = [tcp_x; tcp_y; -tcp_r];

else
    % initialization of the principal point at the center of the image.
    cx = cxcy(1);
    cy = cxcy(2);

    W12 = cx*cy;
    W13 = cx;
    W23 = cy;

    % coefficient matrix
    V = []; % V = [v11,v12,v13,v22,v23,v33];
    b = []; % b = [b11,b12,b13,b22,b23,b33];

    % Eq.(15) in the paper
    vij = @(i,j,iH) [ iH(i,1)*iH(j,1) iH(i,2)*iH(j,2)];
    bij = @(i,j,iH)[iH(j,1)*(W12*iH(i,2)+W13*iH(i,3))+ ...
        iH(j,2)*(W12*iH(i,1)+W23*iH(i,3))+ ...
        iH(j,3)*(iH(i,3)+W13*iH(i,1)+W23*iH(i,2))];

    for ii = 1:nImg

        invH = invHomoSet(:,:,ii);

        v11 = vij(1, 1, invH);
        v12 = vij(1, 2, invH);
        v13 = vij(1, 3, invH);
        v22 = vij(2, 2, invH);
        v23 = vij(2, 3, invH);
        v33 = vij(3, 3, invH);

        b11 = bij(1, 1, invH);
        b12 = bij(1, 2, invH);
        b13 = bij(1, 3, invH);
        b22 = bij(2, 2, invH);
        b23 = bij(2, 3, invH);
        b33 = bij(3, 3, invH);

        V_ii = [[v11;v12;v13;v22;v23;v33],-1/(lambda_i1(ii)^2) * eye(6)];
        b_ii = -[b11;b12;b13;b22;b23;b33];

        V = [V;V_ii];
        b = [b;b_ii];
    end

    wa = V\b;
    w = wa(1:2);
    a = wa(3:8);

    %  Image of the Absolute Conic
    W = [w(1), W12, W13;
        W12, w(2), W23;
        W13, W23, 1];
    cx = W(1,3);
    cy = W(2,3);

    fy = sqrt(W(2,2) - cy*cy);
    fx = sqrt(W(1,1) - cx*cx);

    % initial guess of intrinsic
    K = [fx, 0, cx;
        0, fy, cy;
        0, 0, 1];

    % The position of the camera in the calibration target coordinate system

    A = [a(1), a(2), a(3);
        a(2), a(4), a(5);
        a(3), a(5), a(6)];
    A = A ./ A(3,3);

    % tcp
    tcp_x = A(1,3);
    tcp_y = A(2,3);
    tcp_r = sqrt(A(1,1) - A(1,3)^2);
    tcp = [tcp_x; tcp_y; -tcp_r];
end
