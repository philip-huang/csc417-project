function [barout, COMout, ax] = getBar(COM,w,h,j,d,p,th)
    % Returns the parameters needed to visualize a 2D rectangular bar. 
    % R: orientation of bar
    % p: position of COM of bar 
    % w: width
    % h: height 
    x1 = [-w/2 -w/2 w/2 w/2];
    y1 = [h/2 -h/2 -h/2 h/2];
    x2 = [-j/2 -j/2 j/2 j/2];
    y2 = [-d+j/2 -d-j/2 -d-j/2 -d+j/2];
    x3 = [-j/2 -j/2 j/2 j/2];
    y3 = [d+j/2 d-j/2 d-j/2 d+j/2];
    
    polyin = polyshape({x1,x2,x3},{y1,y2,y3});
    polyout = translate(polyin,p);
    barout = rotate(polyout,rad2deg(th),COM+p); % rotation is always about COM of rigid body
    COMout = COM+p;

    % Define translation frame of each rigid body
    xcol = [1;0];
    ycol = [0;1];
    x = [xcol ycol];
    R = @(theta) [cos(-theta) -sin(-theta);
        sin(-theta) cos(-theta)];
    xcol = R(th)*xcol;
    ycol = R(th)*ycol;
    ax = [xcol ycol];
end