#LANGLADE Maxime
#26/02/20
import os
from math import floor
from os.path import join
import argparse
import subprocess
import cv2
import random

from tqdm import tqdm

#global variables
frame_per_path = 20
nb_group = 8

def extract_frames_in_grp(data_path, output_path):

	os.makedirs(output_path, exist_ok=True)
	for i in range(0,nb_group):
		os.makedirs(join(output_path, 'patch' + str(i)), exist_ok=True)

	reader = cv2.VideoCapture(data_path)
	frame_num = 0
	while (reader.isOpened() | frame_num < nb_group*frame_per_path):
		success, image = reader.read()
		if not success:
			print("frame error for: ", data_path, ' frame n° ', frame_num)
			break

		cv2.imwrite(join(output_path , 'patch' + str(int(frame_num / 20)), '{:04d}.tiff'.format(frame_num)), image)
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
		

	return video


def extract_videos(data_path):

	videos_path = join(data_path, 'videos')
	images_path = join(data_path, 'images')

	sub_path = join(images_path, 'byGroup')
	os.makedirs(sub_path, exist_ok=True)

	nb_video_in_folder = len(os.listdir(videos_path))
	nb_video = 8

	if nb_video > nb_video_in_folder:
		print("Error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	print("nb video=: ", nb_video_in_folder)

	video_ids = [1]


	
	for id in tqdm(video_ids):
		#video = get_video(videos_path , id)
		video = "demo10.mov"
		image_folder = video.split('.')[0]
		extract_frames_in_grp(join(videos_path, video), join(sub_path, image_folder))





if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	args = p.parse_args()

	extract_videos(**vars(args))


