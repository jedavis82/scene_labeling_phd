function [H] = get_hof_masks(py_args)
    close all; 
    addpath(genpath('../hof_matlab_source/'));
    arg_mask_path = py_args(1).arg_mask_path; 
    ref_mask_path = py_args(1).ref_mask_path;
    
    % Load the mask images from disk 
    arg_mask = imread(arg_mask_path);
    ref_mask = imread(ref_mask_path);
    
    % Calculate HOF 
    H(1, :) = hof_raster(arg_mask, ref_mask, 0, 'NumberDirections', 360);
    H(2, :) = hof_raster(arg_mask, ref_mask, 2, 'NumberDirections', 360); 
    H(3, :) = hof_raster(arg_mask, ref_mask, 'hybrid', 'NumberDirections', 360);
end  
