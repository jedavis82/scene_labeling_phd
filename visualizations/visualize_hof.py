"""
For the paper, visualize the histogram results for an image.
"""
import pandas as pd
import cv2
import numpy as np
import matlab.engine
import json
import os
import matplotlib.pyplot as plt
from level_one_utils.hof_giou import compute_hof_display

# Sampling a few image paths for the paper to visualize HOF results with.
ROOT_DIR = '../input/coco_images/'
IMAGE_PATHS = ['000000000036.jpg', '000000000086.jpg', '000000000459.jpg']
HIST_OUTPUT_DIR = './hist_vis/'
if not os.path.exists(HIST_OUTPUT_DIR):
    os.makedirs(HIST_OUTPUT_DIR)

PERSON_L2_SUMMARIES_FILE = '../output/person_level_two_summaries.csv'


def main():
    # Start a matlab instance
    eng = matlab.engine.start_matlab()
    eng.cd('../level_one_utils')

    person_df = pd.read_csv(PERSON_L2_SUMMARIES_FILE, encoding='utf-8', engine='python')

    person_df_filtered = person_df[person_df['relative_path'].isin(IMAGE_PATHS)]

    for idx, row in person_df_filtered.iterrows():
        rel_path = row['relative_path']
        arg_label = row['arg_label']
        arg_box = json.loads(row['arg_bounding_box'])
        ref_label = row['ref_label']
        ref_box = json.loads(row['ref_bounding_box'])
        img_path = ROOT_DIR + rel_path
        orig_img = cv2.imread(img_path, cv2.IMREAD_ANYCOLOR)
        img_width = orig_img.shape[1]
        img_height = orig_img.shape[0]
        matlab_args = {'img_width': img_width,
                       'img_height': img_height,
                       'ref_bounding_box': ref_box,
                       'arg_bounding_box': arg_box
                       }
        f0, f2, hyb = compute_hof_display(eng, matlab_args)

        # Plot the resultant histograms
        hist_file = HIST_OUTPUT_DIR + arg_label + '_' + ref_label
        plt.bar(np.arange(361), f0)
        plt.xticks(np.arange(0, 361, 45))
        plt.savefig(hist_file + '_f0.png')
        plt.show()

        plt.bar(np.arange(361), f2)
        plt.xticks(np.arange(0, 361, 45))
        plt.savefig(hist_file + '_f2.png')
        plt.show()

        plt.bar(np.arange(361), hyb)
        plt.xticks(np.arange(0, 361, 45))
        plt.savefig(hist_file + '_hyb.png')
        plt.show()

        print('Finished creating histogram. Drawing Mask Image')

        blank_img = np.zeros((img_height, img_width, 3), np.uint8)
        cv2.rectangle(blank_img, (arg_box[0], arg_box[1]), (arg_box[2], arg_box[3]), (255, 255, 255), -1)
        cv2.rectangle(blank_img, (ref_box[0], ref_box[1]), (ref_box[2], ref_box[3]), (128, 128, 128), -1)
        img_file = hist_file + '_bounding_boxes.jpg'
        cv2.imwrite(img_file, blank_img)


if __name__ == '__main__':
    main()
