function [K, tcp] = init_intrinsic_param_two_view(points3dSet, points2dSet, nImg, cxcy)

% Initialize the homographies

homoSet = zeros(3,3,nImg);

for ii = 1:nImg

    pts2d =  points2dSet{ii};
    pts3d =  points3dSet{ii};

    H = compute_homography(pts2d, pts3d(1:2,:));
    homoSet(:,:,ii) = H;
end

vij = @(i,j,H) [ H(1,i)*H(1,j);
    H(1,i)*H(2,j) + H(1,j)*H(2,i);
    H(1,i)*H(3,j) + H(1,j)*H(3,i);
    H(2,i)*H(2,j);
    H(2,i)*H(3,j) + H(2,j)*H(3,i);
    H(3,i)*H(3,j)]';


for ii = 1:nImg
    H = homoSet(:,:,ii);
    v11(ii,:) = vij(1, 1, H);
    v12(ii,:) = vij(1, 2, H);
    v13(ii,:) = vij(1, 3, H);
    v22(ii,:) = vij(2, 2, H);
    v23(ii,:) = vij(2, 3, H);
    v33(ii,:) = vij(3, 3, H);
end

% Eq.(25) in the paper

C1 = [v12(1,:); v11(1,:) - v22(1,:); v13(1,:)+v23(1,:)+v33(1,:);
      v12(2,:); v11(2,:) - v22(2,:); v13(2,:)+v23(2,:)+v33(2,:);];

C2 = [zeros(1,6);zeros(1,6);v11(1,:);
      zeros(1,6);zeros(1,6);v11(2,:)];

[V,D]=eig(C1,C2);

ind = find(abs(real(diag(D))) < Inf & not(imag(diag(D))));
[~,ind] = min(abs(diag(D)));

w = V(:,ind);
w = w./w(6);

% initialization of the principal point at the center of the image.
A = [v11*w,zeros(2,2);[0;0],v11*w,[0;0];zeros(2,2),v11*w];
b = [-v13*w;-v23*w;v33*w];

xytcp = A\b;
tcp_x = xytcp(1);
tcp_y = xytcp(2);

if (xytcp(3) - tcp_x^2 - tcp_y^2 < 0)
    error('Error solution, select the other two images.');
end

tcp_r = sqrt(xytcp(3) - tcp_x^2 - tcp_y^2);

tcp = [tcp_x; tcp_y; -tcp_r];

%  Image of the Absolute Conic
W = [w(1), w(2), w(3);
    w(2), w(4), w(5);
    w(3), w(5), w(6)];

cy = (W(1,2)*W(1,3) - W(1,1)*W(2,3)) / (W(1,1)*W(2,2)-W(1,2)*W(1,2));
lambda = W(3,3) - (W(1,3)*W(1,3) + cy*(W(1,2)*W(1,3)-W(1,1)*W(2,3))) / W(1,1);
fx = sqrt(lambda / W(1,1));
fy = sqrt(lambda * W(1,1) / (W(1,1)*W(2,2)-W(1,2)*W(1,2)));
s = -W(1,2) * fx * fx * fy / lambda;
% s = 0;
cx = s * cy / fy - W(1,3) * fx * fx / lambda;

K = [fx, s, cx;
    0, fy, cy;
    0, 0, 1];

end
