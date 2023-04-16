function [r_new, z_new, phi_new] = ray_tracing_cylindrical(r, z, phi, grad_r, grad_z, grad_phi, step_r, step_z, step_phi)
    % Inputs:
    % r, z, phi: current position of ray tracing in cylindrical coordinates
    % grad_r, grad_z, grad_phi: gradient of the function at current position
    % step_size: step size for updating the position of the ray tracing
    % Outputs:
    % r_new, z_new, phi_new: updated position of the ray tracing
    
    % Convert gradient from cylindrical to Cartesian coordinates
    grad_x = grad_r*cos(phi) - grad_phi*sin(phi);
    grad_y = grad_r*sin(phi) + grad_phi*cos(phi);
    
    grad = [grad_x, grad_y, grad_z];

    % Calculate new position of ray tracing in Cartesian coordinates
    [step_x ,  step_y] = pol2cart(step_phi,step_r);
    x_new = r*cos(phi) - grad(1)*step_x;
    y_new = r*sin(phi) - grad(2)*step_y;
    z_new = z - grad(3)*step_z;
%     
    % Convert new position from Cartesian to cylindrical coordinates
    r_new = sqrt(x_new^2 + y_new^2);
    phi_new = atan2(y_new, x_new);
    if (phi_new < -pi) || (phi_new > pi)
        phi_new = wrapToPi(phi_new);
    end
    
    % Return new position of ray tracing in cylindrical coordinates
end