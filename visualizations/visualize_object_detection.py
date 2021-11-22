"""
Visualize the object detection results stored in the CSV file
"""
import pandas as pd
import cv2
import json
import numpy as np

INPUT_IMAGES_DIR = '../input/coco_images/'
OBJECT_DETECTION_FILE = '../input/object_detection.csv'
METADATA_FILE = '../input/metadata.csv'

COLORS = [
    [255, 0, 0],
    [0, 255, 0],
    [0, 0, 255],
    [255, 255, 255]
]


def convert_boxes(json_boxes):
    ret_boxes = []
    boxes = np.array(json.loads(json_boxes))
    for b in boxes:
        x = b[0]
        y = b[1]
        right = b[2]
        bottom = b[3]
        ret_boxes.append([x, y, right, bottom])
    return ret_boxes


def main():
    od_df = pd.read_csv(OBJECT_DETECTION_FILE, encoding='utf-8', engine='python')
    meta_df = pd.read_csv(METADATA_FILE, encoding='utf-8', engine='python')
    for idx, row in od_df.iterrows():
        rel_path = row['relative_path']
        img_meta = meta_df.loc[meta_df['relative_path'] == rel_path]
        img_meta = json.loads(img_meta['labels'].iloc[0])
        img_path = INPUT_IMAGES_DIR + rel_path
        orig_img = cv2.imread(img_path, cv2.IMREAD_COLOR)

        boxes = convert_boxes(row['bounding_boxes'])
        num_objects = row['num_objects']
        if num_objects < 2:
            continue  # No need to show images with 1 detection
        labels = json.loads(row['labels'])
        i_color = 0
        for box, label in zip(boxes, labels):
            color = COLORS[i_color]
            i_color += 1
            if i_color > 3:
                i_color = 0
            cv2.rectangle(orig_img, (box[0], box[1]), (box[2], box[3]), color, 2)
            cv2.putText(orig_img, label, (box[0], box[1]), cv2.FONT_HERSHEY_SIMPLEX, 0.75, color, 2)
        print(f'Image metadata: {img_meta}')
        cv2.imshow('Img', orig_img)
        cv2.waitKey(0)


if __name__ == '__main__':
    main()
