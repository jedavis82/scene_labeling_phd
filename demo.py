"""
Demo script for the scene understanding S2T system
"""

from configparser import ConfigParser
import cv2
import pandas as pd
from object_detection import batch_yolo_detection, draw_detection
from level_one_utils.filter_images import filter_images
from metadata import compute_metadata
from level_one_summaries import compute_level_one_summaries
from level_two_summaries import compute_general_summaries, compute_person_summaries


if __name__ == '__main__':
    config = ConfigParser()
    config.read(['demo_config.ini'])
    test = config['DEFAULT'].get('test')
    model_file = config['DEFAULT'].get('model_file')
    weights_file = config['DEFAULT'].get('weights_file')
    image_dir = config['DEFAULT'].get('image_dir')
    animate_objects_file = config['DEFAULT'].get('animate_objects_file')
    object_ontology_file = config['DEFAULT'].get('ontology_file')
    level_one_output_file = config['DEFAULT'].get('level_one_output_file')

    animate_df = pd.read_csv(animate_objects_file, encoding='utf-8', engine='python')
    ontology_df = pd.read_csv(object_ontology_file, encoding='utf-8', engine='python')

    od_df = batch_yolo_detection(image_dir)
    person_df = filter_images(od_df)
    metadata_df = compute_metadata(person_df, image_dir)
    level_one_df = compute_level_one_summaries(person_df, animate_df)
    general_level_two_df = compute_general_summaries(level_one_df)
    person_level_two_df = compute_person_summaries(level_one_df, ontology_df, metadata_df)
    # Display the final results
    # key, relative_path, img_name, arg_label, arg_bounding_box, ref_label, ref_bounding_box
    # proximity, overlap, f0, f2, hybrid, level_one_summary, person_level_two_summary
    rel_paths = person_level_two_df['relative_path'].unique()
    for p in rel_paths:
        df = person_level_two_df.loc[person_level_two_df['relative_path'] == p]
        img_path = image_dir + p
        img = cv2.imread(img_path, cv2.IMREAD_ANYCOLOR)
        drawn_boxes = []
        for idx, row in df.iterrows():
            arg_box = row['arg_bounding_box']
            ref_box = row['ref_bounding_box']
            arg_label = row['arg_label']
            ref_label = row['ref_label']
            l1 = row['level_one_summary']
            l2 = row['person_level_two_summary']
            print(f'Level one summary: {l1}')
            print(f'Level two summary: {l2}')
            if arg_label not in drawn_boxes:
                cv2.rectangle(img, (arg_box[0], arg_box[1]), (arg_box[2], arg_box[3]), (255, 0, 0), 2)
                cv2.rectangle(img, (ref_box[0], ref_box[1]), (ref_box[2], ref_box[3]), (0, 255, 0), 2)
            if ref_label not in drawn_boxes:
                cv2.putText(img, arg_label, (arg_box[0], arg_box[1]), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (255, 0, 0), 2)
                cv2.putText(img, ref_label, (ref_box[0], ref_box[1]), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (0, 255, 0), 2)
        cv2.imshow('Results', img)
        cv2.waitKey(0)

