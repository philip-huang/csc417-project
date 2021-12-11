function [q0] = generateInitCoords(root, relth, d)
    % for SERIAL arm structures (for now). 
    % root is the world coordinates of the root rigid body element
    % relth is a vector of angular displacements between body elements
    % d is the distance from the COM of the current body to the joint
    % (aligned with the body frame y axis)

    R = @(theta) [cos(theta) -sin(theta);
    sin(theta) cos(theta)];

    n = length(relth);
    
    q0 = root; % first 3 elements is the root

    for i = 1:n    
        q0((i-1)*3+4) = q0((i-1)*3+1) + relth(i);
        q0((i-1)*3+5:(i-1)*3+6)= q0((i-1)*3+2:(i-1)*3+3)  + R(q0((i-1)*3+1))*[0;d] + R(relth(i))*R(q0((i-1)*3+1))*[0;d];
    end 
end