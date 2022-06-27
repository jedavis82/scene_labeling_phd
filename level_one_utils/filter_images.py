"""
Load the object detection CSV file and remove any images that do not contain at least one person detection.
Additionally, remove any images that contain only person detections. They are not useful for level two summaries
"""
import pandas as pd
import json


def remove_counts(labels):
    ret_labels = []
    for l in labels:
        split = l.rsplit('_', 1)[0]
        if split != 'person':
            ret_labels.append(split)
    return ret_labels


def filter_images(df, output_file=None):
    # Store the filtered results as a list of dictionaries and output to a data frame
    results = []
    for idx, row in df.iterrows():
        labels = json.loads(row['labels'])
        if 'person_1' in labels: # Must have person_1 to have any person detections
            no_count_labels = remove_counts(labels)
            if len(no_count_labels) != 0:
                results.append({'relative_path': row['relative_path'], 'bounding_boxes': row['bounding_boxes'],
                                'num_objects': row['num_objects'], 'labels': row['labels'],
                                'img_width': row['img_width'], 'img_height': row['img_height']})
    res_df = pd.DataFrame(results)
    if output_file is not None:
        # Only save to disk if the output file is specified
        res_df.to_csv(output_file, encoding='utf-8', header=True, index=False)
    return res_df
