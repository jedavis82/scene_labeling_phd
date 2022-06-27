"""
Perform object detection using the YOLOv3 model.
    The model is available here: https://pjreddie.com/darknet/yolo/
Send back the results as a pandas data frame
    Each object is numbered according to how many of the same objects are detected in the scene
For this experiment, only images with > 2 object detection results are retained
"""
import cv2
import numpy as np
import os
import pandas as pd
import json
from tqdm import tqdm

# The YOLOv3 configuration model file
MODEL_FILE = './input/models/yolo/yolov3.cfg'
# The YOLOv3 weights file to load into the model
WEIGHTS_FILE = './input/models/yolo/yolov3.weights'
# The list of labels for the coco data set. These are the classifiable objects when using the YOLOv3 model
LABELS_FILE = './input/labels/coco.names'

BOX_CONFIDENCE_THRESHOLD = 0.79  # The cutoff for keeping a bounding box
NMS_THRESHOLD = 0.4  # The cutoff for removing duplicate boxes during non-maxima suppression
INPUT_WIDTH = 416  # YOLOv3 requires a fixed input width and height to the model
INPUT_HEIGHT = 416


def load_label_map():
    labels_dict = {}
    with open(LABELS_FILE, 'r') as lmf:
        i = 0
        for line in lmf:
            labels_dict[i] = line.replace(' ', '_').rstrip('\n')
            i += 1
    return labels_dict


def load_image_paths(images_dir=None):
    image_names = []
    image_paths = []
    for root, dirs, files in os.walk(images_dir):
        for f in files:
            img_path = images_dir + f
            image_paths.append(img_path)
            image_names.append(f)
    # Load the images into an array
    images = []
    for i in image_paths:
        img = cv2.imread(i, cv2.IMREAD_COLOR)
        images.append(img)
    return image_names, images


def get_output_names(net):
    layer_names = net.getLayerNames()
    # return [layer_names[i[0] - 1] for i in net.getUnconnectedOutLayers()]  # Old Python 3.7 and older CV2
    return [layer_names[i - 1] for i in net.getUnconnectedOutLayers()]


def draw_detection(img, boxes, labels):
    for box, label in zip(boxes, labels):
        cv2.rectangle(img, (box[0], box[1]), (box[2], box[3]), (255, 0, 0), 2)
        cv2.putText(img, label, (box[0], box[1]), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (255, 0, 0), 2)
    cv2.imshow("Bounding Box Output", img)
    cv2.waitKey(0)


def batch_yolo_detection(images_dir=None, output_file=None, draw_detections=None):
    assert images_dir is not None, "Must supply input image directory"
    # assert output_file is not None, "Must supply output file path"

    if draw_detections is None:
        draw_detections = False
    labels_dict = load_label_map()
    image_names, images = load_image_paths(images_dir)
    # Load the yolo model using opencv
    net = cv2.dnn.readNetFromDarknet(MODEL_FILE, WEIGHTS_FILE)
    net.setPreferableBackend(cv2.dnn.DNN_BACKEND_OPENCV)
    net.setPreferableTarget(cv2.dnn.DNN_TARGET_CPU)

    # Store each images result as a dict and add to the img_results list
    img_results = []
    for i, n in tqdm(zip(images, image_names), total=len(image_names)):
        # print(f'Processing: {n}')
        blob = cv2.dnn.blobFromImage(i, 1/255, (INPUT_WIDTH, INPUT_HEIGHT), [0, 0, 0], 1, crop=False)
        net.setInput(blob)
        outs = net.forward(get_output_names(net))

        img_height = i.shape[0]
        img_width = i.shape[1]

        class_ids = []
        confidences = []
        boxes = []

        df_boxes = []
        df_labels = []
        df_confidences = []

        # Scan through all results and keep only the ones with the highest confidence scores.
        # Assign the box's class label as the class with the highest score
        for out in outs:
            for detection in out:
                scores = detection[5:]
                class_id = np.argmax(scores)
                confidence = scores[class_id]
                if confidence > BOX_CONFIDENCE_THRESHOLD:
                    center_x = int(detection[0] * img_width)
                    center_y = int(detection[1] * img_height)
                    width = int(detection[2] * img_width)
                    height = int(detection[3] * img_height)
                    left = int(center_x - width / 2)
                    top = int(center_y - height / 2)
                    class_ids.append(class_id)
                    confidences.append(float(confidence))
                    x = left
                    y = top
                    right = (left + width)
                    bottom = (top + height)
                    # Prevent negative coords
                    x = max(0, x)
                    y = max(0, y)
                    right = max(0, right)
                    bottom = max(0, bottom)
                    boxes.append([x, y, right, bottom])

        # Perform NMS to eliminate redundant overlapping bounding boxes with lower confidences
        indices = cv2.dnn.NMSBoxes(boxes, confidences, BOX_CONFIDENCE_THRESHOLD, NMS_THRESHOLD)
        obj_labels = {}
        for idx in indices:
            # i = idx[0]  # Old Python 3.7 and older CV2
            i = idx
            box = boxes[i]
            conf = confidences[i]
            label = labels_dict[int(class_ids[i])]
            # Dedupe the labels
            if label not in obj_labels:
                obj_labels[label] = 1
            else:
                obj_labels[label] += 1
            obj_label = label + '_' + str(obj_labels[label])
            df_labels.append(obj_label)
            df_boxes.append(box)
            df_confidences.append(conf)
        num_objects = len(indices)
        img_results.append({'relative_path': n, 'bounding_boxes': json.dumps(df_boxes), 'num_objects': num_objects,
                            'confidences': df_confidences, 'labels': json.dumps(df_labels),
                            'img_width': img_width, 'img_height': img_height})
        if draw_detections:
            img_path = images_dir + n
            clone = cv2.imread(img_path, cv2.IMREAD_COLOR)
            draw_detection(clone, df_boxes, df_labels)
    # Drop any rows where there are < 2 detections. These can't be processed by the S2T system
    results_df = pd.DataFrame(img_results)
    results_df.drop(results_df[results_df['num_objects'] < 2].index, inplace=True)
    if output_file is not None:
        # Only write to disk when the output file is specified
        results_df.to_csv(output_file, encoding='utf-8', header=True, index=False)
    print('Finished writing results to csv file')
    return results_df
