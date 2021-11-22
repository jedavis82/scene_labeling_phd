function [H] = get_hof(py_args)
    close all; 
    addpath(genpath('./hof_matlab_source/'));
    img_width = py_args(1).img_width + 1; 
    img_height = py_args(1).img_height + 1; 
    arg_bounding_box = py_args(1).arg_bounding_box; 
    ref_bounding_box = py_args(1).ref_bounding_box; 
    
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
    
    % Calculate HOF
    H(1, :) = hof_raster(arg_mask, ref_mask, 0, 'NumberDirections', 360); 
    H(2, :) = hof_raster(arg_mask, ref_mask, 2, 'NumberDirections', 360); 
    H(3, :) = hof_raster(arg_mask, ref_mask, 'hybrid', 'NumberDirections', 360);
end


    