"""
For the paper, generate results for hof for the input images and store in CSV files.
"""
import pandas as pd
import cv2
import numpy as np
import json
import os
from level_one_utils.hof_giou import compute_hof_display

ROOT_DIR = '../input/coco_images/'
IMAGE_PATHS = ['000000000036.jpg', '000000000086.jpg', '000000000459.jpg', '000000000634.jpg', '000000000962.jpg',
               '000000001655.jpg']
HIST_OUTPUT_DIR = './hist_vis/'
if not os.path.exists(HIST_OUTPUT_DIR):
    os.makedirs(HIST_OUTPUT_DIR)

PERSON_L2_SUMMARIES_FILE = '../output/person_level_two_summaries.csv'


def main():
    person_df = pd.read_csv(PERSON_L2_SUMMARIES_FILE, encoding='utf-8', engine='python')

    person_df_filtered = person_df[person_df['relative_path'].isin(IMAGE_PATHS)]

    for idx, row in person_df_filtered.iterrows():
        if row['person_level_two_summary'] == 'Negative Interaction':
            continue  # No need to compute the negative interactions
        rel_path = row['relative_path']
        arg_label = row['arg_label']
        arg_box = json.loads(row['arg_bounding_box'])
        ref_label = row['ref_label']
        ref_box = json.loads(row['ref_bounding_box'])
        img_path = ROOT_DIR + rel_path
        orig_img = cv2.imread(img_path, cv2.IMREAD_ANYCOLOR)
        img_width = orig_img.shape[1]
        img_height = orig_img.shape[0]

        arg_hof = np.zeros((img_width, img_height), dtype='uint8')
        ref_hof = np.zeros((img_width, img_height), dtype='uint8')
        arg_hof[arg_box[0]:arg_box[2], arg_box[1]:arg_box[3]] = 255
        ref_hof[ref_box[0]:ref_box[2], ref_box[1]:ref_box[3]] = 255
        f0, f2, hyb = compute_hof_display(arg_hof, ref_hof)

        out_df = pd.DataFrame(columns=['degrees', 'f0_magnitude', 'f2_magnitude', 'hyb_magnitude'])
        out_df['degrees'] = np.arange(361)

        f0_magnitudes = f0.flatten().tolist()
        f2_magnitudes = f2.flatten().tolist()
        hyb_magnitudes = hyb.flatten().tolist()
        for i in range(361):
            out_df.loc[out_df['degrees'] == i, 'f0_magnitude'] = f0_magnitudes[i]
            out_df.loc[out_df['degrees'] == i, 'f2_magnitude'] = f2_magnitudes[i]
            out_df.loc[out_df['degrees'] == i, 'hyb_magnitude'] = hyb_magnitudes[i]

        out_file = HIST_OUTPUT_DIR + arg_label + '_' + ref_label + '_hists.csv'
        out_df.to_csv(out_file, encoding='utf-8', header=True, index=False)

        # Plot the resultant histograms
        # hist_file = HIST_OUTPUT_DIR + arg_label + '_' + ref_label
        # plt.bar(np.arange(361), f0)
        # plt.xticks(np.arange(0, 361, 45))
        # plt.savefig(hist_file + '_f0.png')
        # plt.show()
        #
        # plt.bar(np.arange(361), f2)
        # plt.xticks(np.arange(0, 361, 45))
        # plt.savefig(hist_file + '_f2.png')
        # plt.show()
        #
        # plt.bar(np.arange(361), hyb)
        # plt.xticks(np.arange(0, 361, 45))
        # plt.savefig(hist_file + '_hyb.png')
        # plt.show()
        #
        # print('Finished creating histogram. Drawing Mask Image')
        #
        # blank_img = np.zeros((img_height, img_width, 3), np.uint8)
        # cv2.rectangle(blank_img, (arg_box[0], arg_box[1]), (arg_box[2], arg_box[3]), (255, 255, 255), -1)
        # cv2.rectangle(blank_img, (ref_box[0], ref_box[1]), (ref_box[2], ref_box[3]), (128, 128, 128), -1)
        # img_file = hist_file + '_bounding_boxes.jpg'
        # cv2.imwrite(img_file, blank_img)


if __name__ == '__main__':
    main()
