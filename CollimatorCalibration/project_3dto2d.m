function pts2D = project_3dto2d(pts3d, Rt, K, d)

if (size(pts3d,1) == 3)
    pts3d(4,:) = ones(1,size(pts3d,2));
end
if (size(pts3d,1) ~=4)
    error('The dimension of wrdpts should be 4 * n. \n');
end

xy = Rt * pts3d;
uv = [K,[0;0;0]] * xy;

xy = xy(1:2,:) ./ [xy(3,:);xy(3,:)];
uv = uv(1:2,:) ./ [uv(3,:);uv(3,:)];

rho = xy(1,:).^2 + xy(2,:).^2;

u0 = K(1,3);
v0 = K(2,3);

% Two-order radial distortion
uv_k(1,:) = uv(1,:) + (uv(1,:) - u0) .* (d(1).*rho + d(2).*rho.^2);
uv_k(2,:) = uv(2,:) + (uv(2,:) - v0) .* (d(1).*rho + d(2).*rho.^2);

pts2D = uv_k(1:2,:);
end