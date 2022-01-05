"""
For the paper, visualize the HOF algorithm output.
Will do these one at a time so that I can place the legends correctly
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HIST_DIR = './hist_vis/'
PLOTTING_FILE = './hist_vis/person_1_umbrella_1_hists.csv'
IMG_OUTPUT_FILE = './hist_vis/umbrella_hof.png'


def main():
    # filepath = Path(HIST_DIR).glob('**/*')
    # hist_files = [x for x in filepath if x.is_file()]

    # The membership function lines to plot
    left_x = [120, 180, 240]
    left_y = [0, 100, 0]
    below_x = [210, 270, 330]
    below_y = [0, 100, 0]
    above_x = [30, 90, 150]
    above_y = [0, 100, 0]
    right1_x = [0, 30, 60]
    right1_y = [0, 100, 0]
    right2_x = [300, 330, 360]
    right2_y = [0, 100, 0]

    df = pd.read_csv(PLOTTING_FILE, encoding='utf-8', engine='python')
    magnitude = list(df['hyb_magnitude'])
    fig, ax = plt.subplots(dpi=300)
    ax.plot(left_x, left_y, color='red', label='left', alpha=0.5)
    plt.fill_between(left_x, left_y, color='red', alpha=0.5)
    ax.plot(below_x, below_y, color='green', label='below', alpha=0.5)
    plt.fill_between(below_x, below_y, color='green', alpha=0.5)
    ax.plot(above_x, above_y, color='blueviolet', label='above', alpha=0.5)
    plt.fill_between(above_x, above_y, color='blueviolet', alpha=0.5)
    ax.plot(right1_x, right1_y, color='darkorange', label='right', alpha=0.5)
    plt.fill_between(right1_x, right1_y, color='darkorange', alpha=0.5)
    ax.plot(right2_x, right2_y, color='darkorange', alpha=0.5)
    plt.fill_between(right2_x, right2_y, color='darkorange', alpha=0.5)

    ax.bar(np.arange(361), magnitude, color='blue')
    ax.set_xticks(np.arange(0, 361, 45))
    plt.legend(loc='upper left')
    plt.xlabel('Degrees')
    plt.ylabel('Magnitude')
    plt.savefig(IMG_OUTPUT_FILE)
    plt.show()


if __name__ == '__main__':
    main()
