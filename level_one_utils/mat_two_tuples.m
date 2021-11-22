function mat_two_tuples(py_args)
    disp(py_args)
    img_width = py_args(1).img_width + 1;
    img_height = py_args(1).img_height + 1; 
    ref_bounding_box = py_args(1).ref_bounding_box; 
    arg_bounding_box = py_args(1).arg_bounding_box;
    ref_label = py_args(1).ref_label;
    arg_label = py_args(1).arg_label;
    
    ref_x = ref_bounding_box{1} + 1;
    ref_y = ref_bounding_box{2} + 1;
    ref_right = ref_bounding_box{3} + 1; 
    ref_bottom = ref_bounding_box{4} + 1;
    
    arg_x = arg_bounding_box{1} + 1; 
    arg_y = arg_bounding_box{2} + 1;
    arg_right =  arg_bounding_box{3} + 1; 
    arg_bottom =  arg_bounding_box{4} + 1; 
    
    ref_mask = zeros(img_width, img_height, 'uint8');
    ref_mask(ref_y:ref_bottom, ref_x:ref_right) = 255;
    arg_mask = zeros(img_width, img_height, 'uint8');
    arg_mask(arg_y:arg_bottom, arg_x:arg_right) = 255; 
    disp(strcat('Reference: ', ref_label))
    disp(strcat('Argument: ', arg_label))
    subplot(1,2,1), imshow(ref_mask)
    subplot(1,2,2), imshow(arg_mask)
end