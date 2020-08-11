#LANGLADE Maxime
#26/02/20
import os
from math import floor
from os.path import join
import argparse
import subprocess
import cv2
import random
import glob

from tqdm import tqdm
from face_extraction import face_extraction

#global variables
frame_per_video = 100
mask_w = 200
mask_h = 200

def extract_frames_training(data_path, output_path, x_center, y_center):

	#print(data_path, " ", x_center, " " ,y_center)

	os.makedirs(output_path, exist_ok=True)

	reader = cv2.VideoCapture(data_path)
	frame_num = 0
	while (reader.isOpened() | frame_num < frame_per_video):
		success, image = reader.read()
		if not success:
			print("frame error for: ", data_path, ' frame n° ', frame_num)
			break

		x_start = x_center - (1/2 * mask_w)
		y_start = y_center - (1/2 * mask_h)
		y_end = int(y_center + mask_h)
		x_end = int(x_center + mask_w)
		patch = image[int(y_start):y_end, int(x_start):x_end]
		cv2.imwrite(join(output_path, '{:04d}.tiff'.format(frame_num)), patch)
		frame_num += 1
	reader.release()		

def get_video(videos_path, id):
	off_set = 0
	success = False

	while not success:
		#print("video id: " + video)

		id_r = (id + off_set) % len(os.listdir(videos_path))
		video = os.listdir(videos_path)[id_r] #todo check unicity in video_ids array

		off_set += 1

		if video[0] == '.':
			continue

		reader = cv2.VideoCapture(join(videos_path, video))
		success, image = reader.read()
		reader.release()

		if not success:
			print("couldn't readvideo n° " + video)
			continue

		x, y = face_extraction(join(videos_path, video)) 

		if x < mask_w or y < mask_h:
			success = False
		

	return video, x, y


def extract_videos(data_path, nb_video):

	videos_path = join(data_path, 'videos')
	images_path = join(data_path, 'images')

	sub_path = join(images_path, 'training')
	os.makedirs(sub_path, exist_ok=True)
	
	'''
	generic_video_name = "id*.mp4"
	video_in_folder = glob.glob(join(video_path, generic_video_name))

	nb_video_in_folder = len(video_in_folder)
	nb_video = int(nb_video)

	if nb_video > nb_video_in_folder:
		print("Error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	video_ids = random.sample(range(1, nb_video_in_folder), nb_video)
'''

	nb_video_in_folder = len(os.listdir(videos_path))
	nb_video = int(nb_video)

	if nb_video > nb_video_in_folder:
		print("Error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	video_ids = random.sample(range(1, nb_video_in_folder), nb_video)

	for id in tqdm(video_ids):
		video, x, y = get_video(videos_path , id)
		image_folder = video.split('.')[0]
		extract_frames_training(join(videos_path, video), join(sub_path, image_folder), x, y)





if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	p.add_argument('--nb_video', '-n', type=str, default='100')
	args = p.parse_args()

	extract_videos(**vars(args))


