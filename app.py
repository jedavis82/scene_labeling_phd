"""
This script serves as the driver file for computing the S2T labels for a set of input images.
Images can be downloaded by using the script in input/import_coco.py
The output of each phase of the S2T system is stored in CSV files for retrieval.
These files are:
    object_detection.csv, person_object_detection.csv, metadata.csv, level_one_summaries.csv,
    general_level_two_summaries.csv, person_level_two_summaries.csv
"""

from object_detection import yolo_object_detection
from level_one_summaries import compute_level_one_summaries
from level_one_utils.filter_images import filter_images
from level_two_summaries import compute_general_summaries, compute_person_summaries
from metadata import compute_metadata
import os
import pandas as pd

IMAGES_DIR = './input/coco_images/'

ANIMATE_OBJECTS_FILE = './input/animate_objects.csv'

OBJECT_ONTOLOGY_FILE = './input/object_ontology.csv'

OBJECT_DETECTION_FILE = './input/object_detection.csv'

PERSON_IMAGES_FILE = './input/person_object_detection.csv'

METADATA_FILE = './input/metadata.csv'

LEVEL_ONE_SUMMARIES_FILE = './output/level_one_summaries.csv'

# Using a boolean flag here because the level one summaries take some time to compute. The script that does this will
# append to the level one summaries file in case it gets interrupted when computing
FINALIZED_LEVEL_ONE = True

GENERAL_LEVEL_TWO_SUMMARIES_FILE = './output/general_level_two_summaries.csv'

PERSON_LEVEL_TWO_SUMMARIES_FILE = './output/person_level_two_summaries.csv'


def main():
    # Compute object localization results
    if not os.path.exists(OBJECT_DETECTION_FILE):
        od_df = yolo_object_detection(IMAGES_DIR, OBJECT_DETECTION_FILE)
    else:
        od_df = pd.read_csv(OBJECT_DETECTION_FILE, encoding='utf-8', engine='python')
    # Can drop the confidence scores from the object detection data frame as they are not used
    od_df.drop(columns=['confidences'], inplace=True)

    # Remove any images from the set of object detection results that does not contain a valid person detection
    if not os.path.exists(PERSON_IMAGES_FILE):
        person_df = filter_images(od_df, PERSON_IMAGES_FILE)
    else:
        person_df = pd.read_csv(PERSON_IMAGES_FILE, encoding='utf-8', engine='python')
    if not os.path.exists(METADATA_FILE):
        metadata_df = compute_metadata(person_df, IMAGES_DIR, METADATA_FILE)
    else:
        metadata_df = pd.read_csv(METADATA_FILE, encoding='utf-8', engine='python')
    # Compute level one summaries for the set of object localization results
    if not FINALIZED_LEVEL_ONE:
        animate_df = pd.read_csv(ANIMATE_OBJECTS_FILE, encoding='utf-8', engine='python')
        animate_objects = list(animate_df['object'])
        l1_df = compute_level_one_summaries(person_df, animate_objects, LEVEL_ONE_SUMMARIES_FILE)
    else:
        l1_df = pd.read_csv(LEVEL_ONE_SUMMARIES_FILE, encoding='utf-8', engine='python')
        print('here')
    # Compute the general domain level two summaries for the level one summaries
    if not os.path.exists(GENERAL_LEVEL_TWO_SUMMARIES_FILE):
        gl2_df = compute_general_summaries(l1_df, GENERAL_LEVEL_TWO_SUMMARIES_FILE)
    else:
        gl2_df = pd.read_csv(GENERAL_LEVEL_TWO_SUMMARIES_FILE, encoding='utf-8', engine='python')
        print('here')
    # Compute the person domain level two summaries for the level one summaries with person - object interactions
    if not os.path.exists(PERSON_LEVEL_TWO_SUMMARIES_FILE):
        ontology_df = pd.read_csv(OBJECT_ONTOLOGY_FILE, encoding='utf-8', engine='python')
        pl2_df = compute_person_summaries(l1_df, ontology_df, metadata_df, PERSON_LEVEL_TWO_SUMMARIES_FILE)
    else:
        pl2_df = pd.read_csv(PERSON_LEVEL_TWO_SUMMARIES_FILE, encoding='utf-8', engine='python')
        print('here')


if __name__ == '__main__':
    main()
