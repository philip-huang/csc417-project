function [allCOM, allBars, allax] = visualizeBar(q,w,h,j,d,m)

    % takes in q (6x1) value, generalized coordinates for each rigid bar
    
    % Rigid-body bar link parameters
    COM_init = [0 0]; % center of mass

    % Define translation frame of each rigid body
    xL = [1;0];
    yL = [0;1];
    x = [xL yL];
    R = @(theta) [cos(theta) -sin(theta);
        sin(theta) cos(theta)];

    numBars = length(q)/3;
    allax = [];
    allBars = [];
    allCOM = [];

    for n = 1:numBars
        p = q(3*(n-1)+2:3*(n-1)+3); % translation
        th = q(3*(n-1)+1); % rotation [rad]
        [bar,COM,ax] =  getBar(COM_init,w,h,j,d,p',th);
        allax = cat(3, allax, ax);
        allBars = [allBars bar];
        allCOM = [allCOM; COM];
    end
    
%     % takes in q (6x1) value, generalized coordinates for each rigid bar
%     
%     % Rigid-body bar link parameters
%     COM_init = [0 0]; % center of mass
% 
%     % Define translation frame of each rigid body
%     xL = [1;0];
%     yL = [0;1];
%     x = [xL yL];
%     R = @(theta) [cos(theta) -sin(theta);
%         sin(theta) cos(theta)];
% 
%     p = q(2:3); % translation
%     th = q(1); % rotation [rad]
%     [bar0,COM0,ax0] =  getBar(COM_init,w,h,j,d,p',th);
%     COM = COM0;
%     %quiver(COM(1), COM(2), ax0(1,1), ax0(1,2),'color',[1 0 0])
%     %quiver(COM(1), COM(2), ax0(2,1), ax0(2,2),'color',[0 1 0])
%     
%     p = q(5:6); % translation
%     th = q(4); % rotation [rad]
%     [bar1,COM1,ax1] =  getBar(COM_init,w,h,j,d,p',th);
%     COM = COM1;
%     %quiver(COM(1), COM(2), ax1(1,1), ax1(1,2),'color',[1 0 0])
%     %quiver(COM(1), COM(2), ax1(2,1), ax1(2,2),'color',[0 1 0])
%    
%     allax = cat(3, ax0,ax1);
%     allBars = [bar0 bar1];
%     allCOM = [COM0; COM1].';
end