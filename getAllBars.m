function [allCOM, allBars, allax] = getAllBars(q,w,h,j,d,m)
    % Takes in q (3*nx1) value, generalized coordinates for each rigid bar
    % where n = number of bars in the whole simulation.

    % Returns the parameters needed to visualize ALL the 2D rectangular bars in the simulation. 
    
    COM_init = [0 0]; % center of mass

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
end