function [] = visualizeAllBars(allCOM, allBars, allax)
    plot(allBars);
    plot(allCOM(:,1),allCOM(:,2),'b o')
    for i = 1:length(allCOM)
        COM = allCOM(i,:);
        ax = allax(:,:,i);
        quiver(COM(1), COM(2), ax(1,1), ax(1,2),'color',[1 0 0])
        quiver(COM(1), COM(2), ax(2,1), ax(2,2),'color',[0 1 0])
    end
    
    numBod = size(allCOM);
    for n = 1:numBod(1)
        labels{n} = uint8(n);
    end 
    text(allCOM(:,1),allCOM(:,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
end