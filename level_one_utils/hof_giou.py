"""
Call the HOF code from Matlab and return the results to the Python script
"""

import numpy as np
from sklearn.preprocessing import minmax_scale
import hofpy


def compute_hof(arg_object=None, ref_object=None, num_directions=360):
    """
    Compute HOF for an object tuple and return the max angles for F0, F2, and Hybrid
    :param arg_object: A binary mask image representing the argument object
    :param ref_object: A binary mask image representing the referrant object
    :param num_directions: The number of histogram of forces directions to compute (default 360)
    :return: The maximum F0, F2, and Hybrid HOF angles
    """
    assert arg_object is not None, "Must supply argument image"
    assert ref_object is not None, "Must supply referrant image"
    f0 = hofpy.F0Hist(arg_object, ref_object, numberDirections=num_directions)
    f2 = hofpy.F2Hist(arg_object, ref_object, numberDirections=num_directions)
    hybrid = hofpy.F02Hist(arg_object, ref_object, numberDirections=num_directions)

    f0 = np.nan_to_num(f0)
    f2 = np.nan_to_num(f2)
    hybrid = np.nan_to_num(hybrid)
    max_f0_angle = np.argmax(f0)
    max_f2_angle = np.argmax(f2)
    max_hybrid_angle = np.argmax(hybrid)
    return max_f0_angle, max_f2_angle, max_hybrid_angle


def compute_hof_display(arg_object=None, ref_object=None, num_directions=360):
    """
    Compute HOF for an object two-tuple and return the histograms for display
    :param arg_object:
    :param ref_object:
    :param num_directions:
    :return:
    """
    assert arg_object is not None, "Must supply the argument mask image"
    assert ref_object is not None, "Must supply the referrant mask image"
    f0 = hofpy.F0Hist(arg_object, ref_object, numberDirections=num_directions)
    f2 = hofpy.F2Hist(arg_object, ref_object, numberDirections=num_directions)
    hybrid = hofpy.F02Hist(arg_object, ref_object, numberDirections=num_directions)
    f0 = np.nan_to_num(f0)
    f2 = np.nan_to_num(f2)
    hybrid = np.nan_to_num(hybrid)
    f0_histograms = minmax_scale(f0, feature_range=(0, 100))
    f2_histograms = minmax_scale(f2, feature_range=(0, 100))
    hyb_histograms = minmax_scale(hybrid, feature_range=(0, 100))

    return f0_histograms.astype(int), f2_histograms.astype(int), hyb_histograms.astype(int)


def compute_giou(arg_obj, ref_obj):
    """
    Compute the generalized intersection over union code described in the paper available from:
    https://giou.stanford.edu/GIoU.pdf Page 4 of the paper has the psuedocode for this algorithm

    :param arg_obj: The argument object represented as a bounding box
    :param ref_obj: The referrant object represented as a bounding box
    :return: The IoU and GIoU scores for the arg and ref objects
    """
    arg_x1 = arg_obj[0]
    arg_y1 = arg_obj[1]
    arg_x2 = arg_obj[2]
    arg_y2 = arg_obj[3]
    area_arg_obj = (arg_x2 - arg_x1) * (arg_y2 - arg_y1)
    # Compute the area of the ref obj
    ref_x1 = ref_obj[0]
    ref_y1 = ref_obj[1]
    ref_x2 = ref_obj[2]
    ref_y2 = ref_obj[3]
    area_ref_obj = (ref_x2 - ref_x1) * (ref_y2 - ref_y1)

    # Calculate the intersection between the arg and ref objects
    x_1_I = max(arg_x1, ref_x1)
    x_2_I = min(arg_x2, ref_x2)
    y_1_I = max(arg_y1, ref_y1)
    y_2_I = min(arg_y2, ref_y2)
    if x_2_I > x_1_I and y_2_I > y_1_I:  # Double check this, I think and is correct here.
        I = (x_2_I - x_1_I) * (y_2_I - y_1_I)
    else:
        I = 0

    # Find the coordinates of the smallest bounding box that encloses both objects
    x_1_c = min(ref_x1, arg_x1)
    x_2_c = max(ref_x2, arg_x2)
    y_1_c = min(ref_y1, arg_y1)
    y_2_c = max(ref_y2, arg_y2)

    # Calculate the area of the smallest enclosing bounding box
    area_b_c = (x_2_c - x_1_c) * (y_2_c - y_1_c)

    # Calculate the IOU (Intersection over Union)
    # IoU = I/U, where U = Area_ref + Area_arg - I
    U = area_arg_obj + area_ref_obj - I
    IoU = I / U

    # Calculate GIoU (Generalized Intersection over Union)
    # GIoU = IoU - (Area_c - U)/Area_c
    GIoU = IoU - ((area_b_c - U) / area_b_c)
    return GIoU, IoU
