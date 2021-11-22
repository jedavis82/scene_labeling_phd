close all; clear all; clc;

load test_data.mat;

%Construct structure array of test images
exRas(10,1) = struct('imgA', [], 'imgB', []);
for testNum = 1:10
    eval(['exRas(testNum).imgA = exRas' num2str(testNum) 'A']);
    eval(['exRas(testNum).imgB = exRas' num2str(testNum) 'B']);
end

for testNum = 1:10
    
    %Show test images
    figure;
    subplot(4,2,1);
    imshow(exRas(testNum).imgA,[]);
    title('Argument');
    subplot(4,2,2);
    imshow(exRas(testNum).imgB,[]);
    title('Referent');
    
    %Calculate HoF
    H(1,:) = hof_raster(exRas(testNum).imgA, exRas(testNum).imgB,0);
    H(2,:) = hof_raster(exRas(testNum).imgA, exRas(testNum).imgB,2);
    H(3,:) = hof_raster(exRas(testNum).imgA, exRas(testNum).imgB,'hybrid');
    
    %Plot histograms
    subplot(4,2,3:4);
    plot(0:2:360,H(1,:));
    axis([0 360 0 max(H(1,:))]);
    title('F0 Histogram');
    subplot(4,2,5:6);
    plot(0:2:360,H(2,:));
    axis([0 360 0 max(H(2,:))]);
    title('F2 Histogram');
    subplot(4,2,7:8);
    plot(0:2:360,H(3,:));
    axis([0 360 0 max(H(3,:))]);
    title('FH Histogram');
    
    
end