#LANGLADE Maxime
#26/02/20
import os
import numpy as np
import dlib
import cv2
from os.path import join
import argparse
import random
from tqdm import tqdm
import glob


from shapely.geometry import Polygon

import json

patch_size = (299, 299)

def rect_to_bb(rect):
	# take a bounding predicted by dlib and convert it
	# to the format (x, y, w, h) as we would normally do
	# with OpenCV
	x = rect.left()
	y = rect.top()
	w = rect.right() - x
	h = rect.bottom() - y
	# return a tuple of (x, y, w, h)
	return (x, y, w, h)


def shape_to_np(shape, dtype="int"):
	# initialize the list of (x, y)-coordinates
	coords = np.zeros((68, 2), dtype=dtype)
	# loop over the 68 facial landmarks and convert them
	# to a 2-tuple of (x, y)-coordinates
	for i in range(0, 68):
		coords[i] = (shape.part(i).x, shape.part(i).y)
	# return the list of (x, y)-coordinates
	return coords

def extract_face(images_folder, nb_frame):
	detector = dlib.get_frontal_face_detector()

	generic_image_name = "f*.png"
	images = glob.glob(join(images_folder, generic_image_name))

	if nb_frame == -1:
		nb_frame = len(images)

	output_path = join(images_folder, "face")
	os.makedirs(output_path, exist_ok=True)

	img_id = 0
	count_img =0

	while count_img < nb_frame:
		img = cv2.imread(images[img_id], flags=cv2.IMREAD_COLOR)
		#if img == None:
		#	print("img is empty")
		#	break

		gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

		# detect faces in the grayscale image
		rects = detector(gray, 1)
		#x1, y1, x2, y2 = face_extraction(frame):

		if len(rects) >= 2:
			print(images_folder)
			print("wrong nb face (", len(rects), ")  in frame ", img_id)
			print("caution try to rect[0]")

		if len(rects) == 0:
			print(images_folder)
			print("NO face (", len(rects), ")  in frame ", img_id)
			print("SKIP")
			img_id = img_id + 1
			continue

		patch = img[rects[0].top():rects[0].bottom(), rects[0].left():rects[0].right()]
		path_resized = cv2.resize(patch, patch_size)
		cv2.imwrite(join(output_path, '{:04d}.tiff'.format(count_img)), path_resized)
		count_img = count_img + 1
		img_id = img_id + 1

	#if len(images) < 20:
	#	print("!!!! ", images_folder, " nb=", len(images))

def process(data_path, nb_frame):
	nb_frame = int(nb_frame)


	generic_video_name = "gen_*"
	images_folders = glob.glob(join(data_path, generic_video_name))

	count = 0
	for id in tqdm(images_folders):
		count = count + 1
		if count >= 24:
			extract_face(id, nb_frame)




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	p.add_argument('--nb_frame', '-n', type=str, default='-1')
	args = p.parse_args()

	process(**vars(args))
