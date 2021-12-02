function [] = visualize(q)
    % takes in q (6x1) value, generalized coordinates for each rigid bar
    
    % Rigid-body bar link parameters
    COM_init = [0 0]; % center of mass
    w = 1; % width of rigid body
    h = 5; % height of rigid body
    j = 0.25; % joint hole size
    d = 2; % distance of the joint holes away from COM, the joint holes on 
           % both sides are the same distance away from the COM
    m = 5; % mass of the object
    
    % Define translation frame of each rigid body
    xL = [1;0];
    yL = [0;1];
    x = [xL yL];
    R = @(theta) [cos(theta) -sin(theta);
        sin(theta) cos(theta)];

p = q(2:3); % translation
th = q(1); % rotation [rad]
[bar0,COM0,ax] =  getBar(COM_init,w,h,j,d,p',th);
COM = COM0;
quiver(COM(1), COM(2), ax(1,1), ax(1,2),'color',[1 0 0])
quiver(COM(1), COM(2), ax(2,1), ax(2,2),'color',[0 1 0])

p = q(5:6); % translation
th = q(4); % rotation [rad]
[bar1,COM1,ax] =  getBar(COM_init,w,h,j,d,p',th);
COM = COM1;
quiver(COM(1), COM(2), ax(1,1), ax(1,2),'color',[1 0 0])
quiver(COM(1), COM(2), ax(2,1), ax(2,2),'color',[0 1 0])

plot([bar0 bar1])
labels = {'A','B'};
allCOM = [COM0; COM1].';
plot(allCOM(1,:),allCOM(2,:),'b o')
text(allCOM(1,:),allCOM(2,:),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

end