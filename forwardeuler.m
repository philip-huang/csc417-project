function [q, qdot] = forwardeuler(q, qdot, qddot, dt)
    % perform next step  
    qdot = qdot + dt*qddot;   
    q = q + dt*qdot;
end