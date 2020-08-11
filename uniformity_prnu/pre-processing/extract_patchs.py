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
frame_per_path = 100
mask_w = 200
mask_h = 200

#functions
def extract_frames_s(data_path, output_path, x1, x2, y):
	os.makedirs(output_path, exist_ok=True)
	os.makedirs(join(output_path, 'patch' + str(0)), exist_ok=True)
	os.makedirs(join(output_path, 'patch' + str(1)), exist_ok=True)
	
	reader = cv2.VideoCapture(data_path)
	frame_num = 0
	while (reader.isOpened() | frame_num < frame_per_path):
		success, image = reader.read()
		if not success:
			break
		patch_1 = image[y:y+mask_h, x1:x1+mask_w]

		patch_2 = image[y:y+mask_h, x2:x2+mask_w]

		cv2.imwrite(join(output_path, 'patch' + str(0), '{:04d}.png'.format(frame_num)),
					 patch_1)

		cv2.imwrite(join(output_path, 'patch' + str(1), '{:04d}.png'.format(frame_num)),
					 patch_2)

		frame_num += 1
	reader.release()

def extract_frames_t(data_path, output_path, x, y):
	os.makedirs(output_path, exist_ok=True)
	os.makedirs(join(output_path, 'patch' + str(0)), exist_ok=True)
	os.makedirs(join(output_path, 'patch' + str(1)), exist_ok=True)

	reader = cv2.VideoCapture(data_path)
	frame_num = 0
	while (reader.isOpened() | frame_num < 2*frame_per_path):
		success, image = reader.read()
		if not success:
			print("frame error for: ", data_path, ' frame n° ', frame_num)
			break

		y_end = int(y + floor(mask_h*1.5))
		x_end = int(x + floor(mask_w*1.5))
		patch = image[int(y):y_end, int(x):x_end]
		if frame_num < frame_per_path:
			cv2.imwrite(join(output_path , 'patch' + str(0), '{:04d}.png'.format(frame_num)),
					 patch)
		else:
			cv2.imwrite(join(output_path , 'patch' + str(1), '{:04d}.png'.format(frame_num - frame_per_path)),
					 patch)
		frame_num += 1
	reader.release()		




def mask_spacial(height, width):
	h_mid = floor(height / 2)
	w_mid = floor(width / 2)

	y = h_mid - floor(mask_h / 2)
	x1 = 15
	x2 = w_mid + 15
	return x1, x2, y


def mask_tamporal(height, width):
	h_mid = floor(height / 2)
	w_mid = floor(width / 2)

	x = w_mid - (3/4 * mask_w)
	y = h_mid - (3/4 * mask_h)

	return x, 0, y
	

def create_mask(video_path, frame ,mode):
	#reader = cv2.VideoCapture(video_path)
	#success, image = reader.read()
	#reader.release()
	#if not success:
	#	print("error during video opening...")
	#	return 0, 0, 0

	height, width, c = frame.shape 

	if mode is 's':
		return mask_spacial(height, width)
	return mask_tamporal(height, width)   	


def get_video(videos_path, id, mode):
	valid_size = False
	off_set = 0
	success = False

	while not success or not valid_size:
		#print("video id: " + video)

		video = os.listdir(videos_path)[id + off_set] #todo check unicity in video_ids array

		off_set += 1

		if video[0] == '.':
			continue

		reader = cv2.VideoCapture(join(videos_path, video))
		success, image = reader.read()
		reader.release()

		if not success:
			print("couldn't readvideo n° " + video)
			continue

		height, width, c = image.shape 
		if height >= (mask_h *1.5) and width > (mask_w *1.5) and not mode is 's':
			#length = int(reader.get(cv2.CAP_PROP_FRAME_COUNT))
			#print("nbframe= ", length)
			#if (length >= frame_per_path * 2):
			valid_size = True
		elif height >= mask_h and width > (2 * mask_w + 30) and mode is 's':	
			valid_size = True	

	return video, image


def extract_videos(data_path, nb_video, mode):

	videos_path = join(data_path, 'videos')
	images_path = join(data_path, 'images')

	if mode is 's':
		sub_path = join(images_path, 'spacial')
		os.makedirs(sub_path, exist_ok=True)
	else:
		sub_path = join(images_path, 'temporal')
		os.makedirs(sub_path, exist_ok=True)

	nb_video_in_folder = len(os.listdir(videos_path))
	nb_video = int(nb_video)

	if nb_video > nb_video_in_folder:
		print("error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	video_ids = random.sample(range(1, nb_video_in_folder), nb_video)
	
	for id in tqdm(video_ids):
		video, frame = get_video(videos_path , id, mode)
		image_folder = video.split('.')[0]
		x1, x2, y = create_mask(join(videos_path, video), frame, mode)    
		
		if mode is 's':
			#SPACIAL
			extract_frames_s(join(videos_path, video), join(sub_path, image_folder), x1, x2, y)
		else:
			#TEMPORAL
			extract_frames_t(join(videos_path, video), join(sub_path, image_folder), x1, y)





if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	p.add_argument('--nb_video', '-n', type=str, default='30')
	p.add_argument('--mode', '-m', type=str)
	args = p.parse_args()

	extract_videos(**vars(args))


