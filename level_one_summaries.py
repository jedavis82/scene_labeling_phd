"""
Generate the level one summaries for the object localization results.
This method uses the Histogram of Forces algorithm to compute the spatial relationship information.
This method uses the Generalized Intersection Over Union method to compute the proximity and overlap information.
Return a pandas data frame containing the level one summaries for each image
    Each level one summary will have the proximity, overlap, and spatial relationship numerical values
    As well as the level one summary labels for each of the object two-tuples

The Matlab engine must be installed on the end user's machine. To do this follow these instructions from your
virtual environment:
https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
"""
import matlab.engine
import pandas as pd
import json
from itertools import permutations
from tqdm import tqdm
import numpy as np
from level_one_utils.hof_giou import compute_hof, compute_giou
from level_one_utils.defuzz_level_one import Defuzz
import os


DIRECTION_LOOKUP = {'Right': 'is to the right of',
                    'Above Right': 'is above and to the right of',
                    'Above': 'is above',
                    'Above Left': 'is above and to the left of',
                    'Left': 'is to the left of',
                    'Below Left': 'is below and to the left of',
                    'Below': 'is below',
                    'Below Right': 'is below and to the right of'}

PROCESSED_IMAGES_FILE = './level_one_utils/processed_images.txt'

L1_RESULTS = []


def convert_boxes(json_boxes):
    ret_boxes = []
    boxes = np.array(json.loads(json_boxes))
    for b in boxes:
        x = b[0]
        y = b[1]
        right = b[2]
        bottom = b[3]
        # Must convert from numpy int32 data type to python's int or the matlab engine will throw an error
        ret_boxes.append([int(x), int(y), int(right), int(bottom)])
    return ret_boxes


def get_two_tuples(img_labels=None):
    """
    Compute the two tuple permutations for the objects in an image.
    Order the tuples by Person -> Animate Object -> Inanimate Object
    :param img_labels: Labels for an input image's object detection results
    :return: Pairs of two-tuples without inverse relationships
    """
    two_tuples = list(permutations(img_labels, 2))
    no_inverse_tuples = []
    no_inverse_key = []
    for tup in two_tuples:
        check = tup[1] + tup[0]
        if check in no_inverse_key:
            continue  # No processing inverse relationships
        no_inverse_tuples.append(tup)
        no_inverse_key.append(tup[0] + tup[1])
    return no_inverse_tuples


def order_arg_ref_pair(al, rl, animate_objects):
    # Use the animate objects list to select the order of objects
    # Person > Animate > Inanimate
    al_split = '_'.join(al.split('_')[:-1])
    rl_split = '_'.join(rl.split('_')[:-1])
    if al_split == 'person':
        return al, rl  # People always take the precedent
    if rl_split == 'person':
        return rl, al  # People always take the precedent
    if al_split in animate_objects:
        return al, rl  # Animate objects take precedent over inanimate objects
    if rl_split in animate_objects:
        return rl, al  # Animate objects take precedent over inanimate objects
    else:
        return al, rl  # Both are inanimate objects


def get_consensus_angle(f0, f2, hyb):
    if f0 == f2:
        return f0
    if f0 == hyb:
        return hyb
    if f2 == hyb:
        return hyb
    else:
        return hyb


def construct_level_one_summary(arg_label, ref_label, overlap_label, sr_label):
    direction = DIRECTION_LOOKUP[sr_label]
    if overlap_label == 'Overlap':
        return f'{arg_label} overlaps and {direction} {ref_label}'
    else:
        return f'{arg_label} {direction} {ref_label}'


def write_to_disk(processed_images, level_one_file):
    with open(PROCESSED_IMAGES_FILE, 'a') as f:
        for p in processed_images:
            f.write(p + '\n')
    if not os.path.exists(level_one_file):
        header = True
        mode = 'w'
    else:
        header = False
        mode = 'a'
    df = pd.DataFrame(L1_RESULTS)
    df.to_csv(level_one_file, mode=mode, header=header, index=False, encoding='utf-8')


def compute_level_one_summaries(od_df, animate_objects, output_file):
    assert od_df is not None, 'Must supply object detection results'
    assert animate_objects is not None, 'Must supply list of animate objects'
    assert output_file is not None, 'Must supply an output file'

    if os.path.exists(PROCESSED_IMAGES_FILE):
        with open(PROCESSED_IMAGES_FILE, 'r') as f:
            processed_images = set(f.read().splitlines())
    else:
        processed_images = set()

    eng = matlab.engine.start_matlab()
    eng.cd('./level_one_utils/')
    defuzzer = Defuzz()
    # od_df = od_df[:50]  # Used for testing purposes

    num_processed = 0
    for idx, row in tqdm(od_df.iterrows(), total=len(od_df)):
        boxes = convert_boxes(row['bounding_boxes'])
        labels = json.loads(row['labels'])
        img_width = row['img_width']
        img_height = row['img_height']
        rel_path = row['relative_path']
        img_name = rel_path.rsplit('.', 1)[0]
        if rel_path in processed_images:
            continue  # Already processed this image
        # Map the bounding boxes to their corresponding labels for use when computing the level one summaries
        label_box_map = {}
        for box, label in zip(boxes, labels):
            label_box_map[label] = box

        img_two_tuples = get_two_tuples(img_labels=labels)
        for tup in img_two_tuples:
            arg_label, ref_label = order_arg_ref_pair(tup[0], tup[1], animate_objects)
            arg_box = label_box_map[arg_label]
            ref_box = label_box_map[ref_label]
            key = f'{img_name}_{arg_label}_{ref_label}'

            giou, iou = compute_giou(arg_box, ref_box)
            matlab_args = {'img_width': img_width,
                           'img_height': img_height,
                           'ref_bounding_box': ref_box,
                           'arg_bounding_box': arg_box}
            f0, f2, hybrid = compute_hof(eng, matlab_args)
            sr_angle = get_consensus_angle(f0, f2, hybrid)
            overlap_label, sr_label = defuzzer.defuzzify_results(iou, sr_angle)
            l1_summary = construct_level_one_summary(arg_label, ref_label, overlap_label, sr_label)
            L1_RESULTS.append({
                'key': key, 'relative_path': rel_path, 'img_name': img_name, 'arg_label': arg_label,
                'arg_bounding_box': arg_box, 'ref_label': ref_label, 'ref_bounding_box': ref_box,
                'overlap': iou, 'proximity': giou, 'f0': f0, 'f2': f2, 'hybrid': hybrid,
                'level_one_summary': l1_summary
            })

            processed_images.add(rel_path)
            num_processed += 1
            if num_processed % 10 == 0:
                write_to_disk(processed_images, output_file)
                processed_images.clear()
                L1_RESULTS.clear()

    # Final write to disk after the loop finishes
    write_to_disk(processed_images, output_file)
    results = pd.read_csv(output_file, encoding='utf-8', engine='python')
    return results

