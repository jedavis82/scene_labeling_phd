function [H] = get_hof_display(py_args)
    close all;
    addpath(genpath('../hof_matlab_source/'));
    img_width = py_args(1).img_width + 1; 
    img_height = py_args(1).img_height + 1; 
    arg_bounding_box = py_args(1).arg_bounding_box; 
    ref_bounding_box = py_args(1).ref_bounding_box; 
    arg_label = py_args(1).arg_label;
    ref_label = py_args(1).ref_label;
    hist_output_dir = py_args(1).hist_output_dir; 
    
    arg_x = arg_bounding_box{1} + 1;
    arg_y = arg_bounding_box{2} + 1; 
    arg_right = arg_bounding_box{3} + 1; 
    arg_bottom = arg_bounding_box{4} + 1; 
    
    if arg_x > img_height
        arg_x = img_height;
    end
    if arg_y > img_width 
        arg_y = img_width;
    end
    if arg_right > img_height
        arg_right = img_height;
    end
    if arg_bottom > img_width
        arg_bottom = img_width;
    end
    
    
    ref_x = ref_bounding_box{1} + 1; 
    ref_y = ref_bounding_box{2} + 1;
    ref_right = ref_bounding_box{3} + 1; 
    ref_bottom = ref_bounding_box{4} + 1; 
    
    if ref_x > img_height
        ref_x = img_height;
    end
    if ref_y > img_width
        ref_y = img_width;
    end
    if ref_right > img_height
        ref_right = img_height;
    end
    if ref_bottom > img_width
        ref_bottom = img_width; 
    end
    
    arg_mask = zeros(img_width, img_height, 'uint8');
    arg_mask(arg_y:arg_bottom, arg_x:arg_right) = 255;
    ref_mask = zeros(img_width, img_height, 'uint8');
    ref_mask(ref_y:ref_bottom, ref_x:ref_right) = 255; 
    
    mask_img = zeros(img_width, img_height, 'uint8');
    mask_img(arg_y:arg_bottom, arg_x:arg_right) = 128;
    mask_img(ref_y:ref_bottom, ref_x:ref_right) = 255;
    
    % disp(size(arg_mask));
    % disp(size(ref_mask));
    
    fig = figure; 
    subplot(4, 2, 1:2); 
    imshow(mask_img, []); 
    title('Arg(gray) Ref(white)');
%     subplot(4, 2, 2); 
%     imshow(ref_mask, []);
%     title(strcat('Referent: ', ref_label));
    
    % Calculate HOF
    H(1, :) = hof_raster(arg_mask, ref_mask, 0, 'NumberDirections', 360); 
    H(2, :) = hof_raster(arg_mask, ref_mask, 2, 'NumberDirections', 360); 
    H(3, :) = hof_raster(arg_mask, ref_mask, 'hybrid', 'NumberDirections', 360);
    
    % Plot histograms 
    subplot(4,2,3:4);
    plot(0:1:360,H(1,:)); 
    axis([0 360 0 max(H(1,:))]);
    title('F0 Histogram');
    subplot(4,2,5:6);
    plot(0:1:360,H(2,:));
    axis([0 360 0 max(H(2,:))]);
    title('F2 Histogram');
    subplot(4,2,7:8);
    plot(0:1:360,H(3,:));
    axis([0 360 0 max(H(3,:))]);
    title('FH Histogram');
    out_file = strcat(hist_output_dir, arg_label, '_', ref_label, '_histograms.png');
    saveas(fig, out_file);
end

% I'm using the 360 directions so I can get an exact 
% measurement for each angle. May slow it down a good bit
% so need to figure out how to use the 180 directions
