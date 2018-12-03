import numpy as np
import cv2


def find_barycenter(thresholded_img):
    """
    Find the barycenter using moments
    :param thresholded_img: segmented image
    :return: 2 elements list with the barycenter coordinates
    """
    contours = cv2.findContours(thresholded_img, 1, 2)
    cnt = contours[0]
    m = cv2.moments(cnt)
    cx = int(m['m10'] / m['m00'])
    cy = int(m['m01'] / m['m00'])
    return [cx, cy]


def find_weighted_centroid(img, mask):
    """
    :param img: numpy ndarray
    :param mask: segmented image with one blob
    :return: 2 elements list with the centroid coordinates
    """
    weighted_blob = np.multiply(img, mask)
    int_sum = np.sum(weighted_blob)

    rows, cols = img.shape

    area = np.count_nonzero(weighted_blob)
    radius = np.sqrt(area/np.pi)

    # sum by rows
    a = np.matrix(range(rows)).transpose()
    b = np.tile(a, (1, rows))
    y = np.sum(np.multiply(weighted_blob, b)) / int_sum

    # sum by columns
    c = np.matrix(range(cols))
    d = np.tile(c, (cols, 1))
    x = np.sum(np.multiply(weighted_blob, d)) / int_sum

    return [int(round(x)), int(round(y))], area, radius, weighted_blob
