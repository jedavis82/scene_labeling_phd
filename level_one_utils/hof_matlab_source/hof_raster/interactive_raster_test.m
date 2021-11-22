close all; clear all; clc;

imgA = zeros(20);
imgB = zeros(20);

numDir = 180;
fuzzLevel = 255;

while(1)
    dispImg = zeros(20,20,3);
    dispImg(:,:,1) = imgA;
    dispImg(:,:,3) = imgB;

    figure(1);
    %imshow(dispImg, 'InitialMagnification', 'fit');
    image(dispImg);
    axis off;
    axis image;
        
    if max(imgA(:)) > 0 && max(imgB(:)) > 0
        H = zeros(3,numDir+1);
        H(1,:) = hof_raster(imgA*fuzzLevel,imgB*fuzzLevel,0,'NumberDirections',numDir);
        H(2,:) = hof_raster(imgA*fuzzLevel,imgB*fuzzLevel,2,'NumberDirections',numDir);
        H(3,:) = hof_raster(imgA*fuzzLevel,imgB*fuzzLevel,'hybrid','NumberDirections',numDir);

        figure(2);
        if(max(isnan(H(1,:))) == 0)
            subplot(3,1,1);
            plot(0:360/numDir:360,H(1,:));
            axis([0,360,0,max(max(H(1,:)), 1)]);
            title('F0 Histogram');
        end
        if(max(isnan(H(2,:))) == 0)
            subplot(3,1,2);
            plot(0:360/numDir:360,H(2,:));
            axis([0,360,0,max(max(H(2,:)), 1)]);
            title('F2 Histogram');
        end
        if(max(isnan(H(3,:))) == 0)
            subplot(3,1,3);
            plot(0:360/numDir:360,H(3,:));
            axis([0,360,0,max(max(H(3,:)),1)]);
            title('F02 Histogram');
        end
    end

    figure(1);
    [x, y, button] = ginput(1);
    if(button == 1)
        if(imgA(round(y),round(x)) == 0)
            imgA(round(y),round(x)) = 1;
        else
            imgA(round(y),round(x)) = 0;
        end
    elseif(button == 3)
        if(imgB(round(y),round(x)) == 0)
            imgB(round(y),round(x)) = 1;
        else
            imgB(round(y),round(x)) = 0;
        end
    end
end