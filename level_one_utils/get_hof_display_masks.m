function [H] = get_hof_display_masks(py_args)
    close all; 
    addpath(genpath('../hof_matlab_source/'));
    arg_mask_path = py_args(1).arg_mask_path; 
    ref_mask_path = py_args(1).ref_mask_path; 
    arg_label = py_args(1).arg_label;
    ref_label = py_args(1).ref_label;
    hist_output_dir = py_args(1).hist_output_dir;
    
    % Load the mask images from disk 
    arg_mask = imread(arg_mask_path); 
    ref_mask = imread(ref_mask_path);
    
    fig = figure; 
    subplot(5, 2, 1:2); 
    imshow(arg_mask, []);
    title('Arg object');
    subplot(5, 2, 3:4);
    imshow(ref_mask, []);
    title('Ref object');
    
    % Calculate HOF 
    H(1, :) = hof_raster(arg_mask, ref_mask, 0, 'NumberDirections', 360); 
    H(2, :) = hof_raster(arg_mask, ref_mask, 2, 'NumberDirections', 360); 
    H(3, :) = hof_raster(arg_mask, ref_mask, 'hybrid', 'NumberDirections', 360);
    
     % Plot histograms 
    subplot(5,2,5:6);
    plot(0:1:360,H(1,:)); 
    axis([0 360 0 max(H(1,:))]);
    title('F0 Histogram');
    subplot(5,2,7:8);
    plot(0:1:360,H(2,:));
    axis([0 360 0 max(H(2,:))]);
    title('F2 Histogram');
    subplot(5,2,9:10);
    plot(0:1:360,H(3,:));
    axis([0 360 0 max(H(3,:))]);
    title('FH Histogram');
    
    out_file = strcat(hist_output_dir, arg_label, '_', ref_label, '_histograms.png');
    saveas(fig, out_file);
end

    
    