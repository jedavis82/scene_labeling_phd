"""
Compute the meta data for each image using the Inception model available from PyTorch
"""
import pandas as pd
import torch
from torchvision import models
import cv2
import numpy as np
import os
import json
from tqdm import tqdm

LABELS_FILE = './input/labels/ilsvrc2012_wordnet_lemmas.txt'
TORCH_HOME_DIR = './input/models/torchvision_models/resnet/'

# Constants used by the resnet model in PyTorch
IMAGE_SIZE = 224
MEAN = [0.485, 0.456, 0.406]
STD = [0.229, 0.224, 0.225]


def preprocess_image(img):
    image = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    image = cv2.resize(image, (IMAGE_SIZE, IMAGE_SIZE))
    image = image.astype("float32") / 255.0
    image -= MEAN
    image /= STD
    image = np.transpose(image, (2, 0, 1))
    image = np.expand_dims(image, 0)
    return image


def compute_metadata(df, images_dir, output_file=None):
    os.environ['TORCH_HOME'] = TORCH_HOME_DIR
    imagenet_labels = dict(enumerate(open(LABELS_FILE)))
    img_paths = list(df['relative_path'].unique())
    model = models.resnet50(pretrained=True)
    model = model.to('cpu')
    model.eval()
    results = []

    for i in tqdm(img_paths, total=len(img_paths)):
        img_path = images_dir + i
        img = cv2.imread(img_path, cv2.IMREAD_COLOR)
        image = preprocess_image(img)
        image = torch.from_numpy(image)
        image = image.to('cpu')
        logits = model(image)
        probabilities = torch.nn.Softmax(dim=1)(logits)
        sorted_probabilities = torch.argsort(probabilities, dim=1, descending=True)

        df_labels = []
        df_confs = []

        for (_, idx) in enumerate(sorted_probabilities[0, :5]):
            _label = imagenet_labels[idx.item()].strip()
            _conf = float(probabilities[0, idx.item()])
            df_labels.append(_label)
            df_confs.append(_conf)
        results.append({
            'relative_path': i,
            'labels': json.dumps(df_labels),
            'confidences': df_confs,
            'num_labels': len(df_labels)
        })

    results_df = pd.DataFrame(results)
    if output_file is not None:
        # Only save to csv if output file is specified
        results_df.to_csv(output_file, encoding='utf-8', header=True, index=False)
    return results_df

