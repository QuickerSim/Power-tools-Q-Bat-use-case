function [grad_u, grad_u_faces] = calculateUnsFVGradients(fv_grid, u, n_fv_gradient_subiters, grad_u)

fvmesh = fv_grid.fvmesh;
grad_u_faces = zeros(3, size(fvmesh.face_con, 2));
vk = zeros(size(grad_u_faces));

internal_faces = fvmesh.face_con(4,:) > 0 & fvmesh.face_con(5,:) > 0;
boundary_faces = ~internal_faces;

tmpA = fvmesh.face_intersect - fvmesh.cell_centroid(:,fvmesh.face_con(4,:));
% dAf = sqrt(tmp(1,:).^2 + tmp(2,:).^2 + tmp(3,:).^2);
dAf = vecnorm(tmpA);
tmpB = zeros(size(tmpA));
tmpB(:, internal_faces) = fvmesh.face_intersect(:, internal_faces) ...
    - fvmesh.cell_centroid(:,fvmesh.face_con(5, internal_faces));
dBf = vecnorm(tmpB);

u_face_intersect = zeros(1, size(fvmesh.face_con, 2));

u_face_intersect(internal_faces) = u(fvmesh.face_con(4, internal_faces)) .* ...
    dAf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces)) + ...
    u(fvmesh.face_con(5, internal_faces)) .* ...
    dBf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces));

u_face_intersect(boundary_faces) = u(fvmesh.face_con(4, boundary_faces)) + ...
    (dot(grad_u(:, fvmesh.face_con(4,boundary_faces)), tmpA(:, boundary_faces)));


for j = 1 : n_fv_gradient_subiters
    
    %skewness correction
    u_face_centroid = u_face_intersect + dot(grad_u_faces, fvmesh.face_centroid - fvmesh.face_intersect);

    grad_u = zeros(3, size(fvmesh.cell_con, 2));

    vk(1, :) = u_face_centroid .* fvmesh.normal(1, :) .* fvmesh.face_area;
    vk(2, :) = u_face_centroid .* fvmesh.normal(2, :) .* fvmesh.face_area;
    vk(3, :) = u_face_centroid .* fvmesh.normal(3, :) .* fvmesh.face_area;

    for i=1:size(fvmesh.face_con, 2)
        grad_u(:, fvmesh.face_con(4, i)) = grad_u(:, fvmesh.face_con(4, i)) + vk(:, i);
        if fvmesh.face_con(5, i) ~= 0
            grad_u(:, fvmesh.face_con(5, i)) = grad_u(:, fvmesh.face_con(5, i)) - vk(:, i);
        end
    end
    grad_u(1, :) = grad_u(1, :) ./ fvmesh.cell_area;
    grad_u(2, :) = grad_u(2, :) ./ fvmesh.cell_area;
    grad_u(3, :) = grad_u(3, :) ./ fvmesh.cell_area;
 
%     grad_u_faces = zeros(3, size(fvmesh.face_con, 2));
    grad_u_faces(1, internal_faces) = grad_u(1, fvmesh.face_con(4, internal_faces)) .* ...
        dAf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces)) + ...
        grad_u(1, fvmesh.face_con(5, internal_faces)) .* ...
        dBf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces));
    
    grad_u_faces(2, internal_faces) = grad_u(2, fvmesh.face_con(4, internal_faces)) .* ...
        dAf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces)) + ...
        grad_u(2, fvmesh.face_con(5, internal_faces)) .* ...
        dBf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces));
    
    grad_u_faces(3, internal_faces) = grad_u(3, fvmesh.face_con(4, internal_faces)) .* ...
        dAf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces)) + ...
        grad_u(3, fvmesh.face_con(5, internal_faces)) .* ...
        dBf(internal_faces) ./ (dAf(internal_faces) + dBf(internal_faces));
    
    grad_u_faces(1, boundary_faces) = grad_u(1, fvmesh.face_con(4, boundary_faces));
    grad_u_faces(2, boundary_faces) = grad_u(2, fvmesh.face_con(4, boundary_faces));
    grad_u_faces(3, boundary_faces) = grad_u(3, fvmesh.face_con(4, boundary_faces));
    
end

end

