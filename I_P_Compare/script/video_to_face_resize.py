#LANGLADE Maxime
#19/06/20
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

output_main = "/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/I_P_Compare/all_resize"

patch_size = (299, 299)

def extract_face(video_path, nb_frame, sub_folder):
	detector = dlib.get_frontal_face_detector()

	video_name = video_path.split("/")[-1]
	video_name = video_name.split(".")[0]

	output_path = join(output_main, sub_folder, video_name)
	os.makedirs(output_path, exist_ok=True)

	reader = cv2.VideoCapture(video_path)

	img_id = 0

	#reader.isOpened()

	while img_id < nb_frame:
		success, img = reader.read()

		if not success:
			print("error reading frame", video_name)
			continue

		gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

		# detect faces in the grayscale image
		rects = detector(gray, 1)
		#x1, y1, x2, y2 = face_extraction(frame):

		if len(rects) >= 2:
			print(video_name)
			print("wrong nb face (", len(rects), ")  in frame ", img_id)
			print("caution try to rect[0]")

		if len(rects) == 0:
			print(video_name)
			print("NO face (", len(rects), ")  in frame ", img_id)
			print("SKIP")
			continue

		patch = img[rects[0].top():rects[0].bottom(), rects[0].left():rects[0].right()]
		path_resized = cv2.resize(patch, patch_size)
		cv2.imwrite(join(output_path, '{:04d}.tiff'.format(img_id)), path_resized)
		img_id = img_id + 1


def get_videos_path(images_folders):

	videos_path_list = []

	for im_folder in images_folders:
		video_name = im_folder.split("/")[-1]
		video_name = video_name + ".mp4"
		if video_name.count("id") == 1:
			video_path = join("/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/videos", video_name)
		elif  video_name.count("id") == 2:
			video_path = join("/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/videos", video_name)
		else:
			print("error count id of :" + video_name + " (" + video.count("id") + ")")
			exit()

		videos_path_list.append(video_path)

	return videos_path_list


def process(data_path):


	generic_video_name = "id*"
	images_folders_test = glob.glob(join(data_path, "test" ,generic_video_name))
	images_folders_train = glob.glob(join(data_path, "train_val" ,generic_video_name))

	video_test_path = get_videos_path(images_folders_test)
	video_train_path = get_videos_path(images_folders_train)

	print("extract test base frames")
	for video in tqdm(video_test_path):
		extract_face(video, 15, "test")

	print("extract train base frames")
	for video in tqdm(video_train_path):
		extract_face(video, 20, "train_val")




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	args = p.parse_args()

	process(**vars(args))
