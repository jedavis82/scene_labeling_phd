"""
Visualize the level one summary results stored in the CSV file
"""
import pandas as pd
import cv2
import json
import numpy as np
import os

INPUT_IMAGES_DIR = '../input/coco_images/'
LEVEL_ONE_SUMMARIES_FILE = '../input/level_one_summaries.csv'

OUTPUT_DIR = '../input/image_results/'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

COLORS = [
    [255, 0, 0],
    [0, 255, 0],
    [0, 0, 255],
    [255, 255, 255]
]


def convert_box(json_box):
    box = np.array(json.loads(json_box))
    return [int(box[0]), int(box[1]), int(box[2]), int(box[3])]


def draw_img_results(img, boxes, labels):
    clone = img.copy()
    i_color = 0
    for b, l in zip(boxes, labels):
        color = COLORS[i_color]
        i_color += 1
        if i_color > 3:
            i_color = 0
        cv2.rectangle(clone, (b[0], b[1]), (b[2], b[3]), color, 2)
        if b[0] - 10 < 0:
            l_x = b[0] + 10
        else:
            l_x = b[0] + 20
        if b[1] - 20 < 0:
            l_y = b[1] + 20
        else:
            l_y = b[1] - 10
        cv2.putText(clone, l, (l_x, l_y), cv2.FONT_HERSHEY_SIMPLEX, 0.75, color, 2)
    return clone


def main():
    l1_df = pd.read_csv(LEVEL_ONE_SUMMARIES_FILE, encoding='utf-8', engine='python')
    img_paths = list(l1_df['relative_path'].unique())
    for p in img_paths:
        img_df = l1_df.loc[l1_df['relative_path'] == p]
        img_path = INPUT_IMAGES_DIR + p
        orig_img = cv2.imread(img_path, cv2.IMREAD_COLOR)
        boxes = []
        labels = []
        l1_summaries = []
        for idx, row in img_df.iterrows():
            arg_label = row['arg_label']
            if not arg_label in labels:
                labels.append(arg_label)
                boxes.append(convert_box(row['arg_bounding_box']))
            ref_label = row['ref_label']
            if not ref_label in labels:
                labels.append(ref_label)
                boxes.append(convert_box(row['ref_bounding_box']))
            l1_summary = row['level_one_summary']
            if not l1_summary in l1_summaries:
                l1_summaries.append(row['level_one_summary'])
        result_img = draw_img_results(orig_img, boxes, labels)
        for l in l1_summaries:
            print(l)
        print()
        cv2.imshow('Level One Summaries', result_img)
        in_key = cv2.waitKey(0)
        # If 's' is pressed, store the image to the output directory
        if in_key == ord('s'):
            output_path = OUTPUT_DIR + p
            cv2.imwrite(output_path, result_img)


if __name__ == '__main__':
    main()
