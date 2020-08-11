#LANGLADE Maxime
#16/04/20
import os
from math import floor
from os.path import join
import argparse
import subprocess
import cv2
import random
import glob

from tqdm import tqdm

def extract_I_frames(video_path, output_path):

	#video_path = "/Volumes/VERBATIM\ HD/" + video_path[21:]
	#print("V_P =", video_path)

	os.makedirs(output_path, exist_ok=True)


	out_file = join(output_path, '%03d.png')
	select_cmd = " -f image2 -vf " + '''"select='eq(pict_type,PICT_TYPE_I)'"'''
	cmd = "ffmpeg -i " + video_path + select_cmd + " -vsync vfr " + out_file
	
	os.system(cmd)              


def process(data_path, nb_video):

	videos_path = join(data_path, 'videos')
	images_path = join(data_path, 'images')

	#sub_path = join(images_path, 'trainingV2')
	sub_path = join(images_path, 'testV2')
	os.makedirs(sub_path, exist_ok=True)


	generic_video_name = "id*.mp4"
	video_in_folder = glob.glob(join(videos_path, generic_video_name))

	nb_video_in_folder = len(video_in_folder)
	nb_video = int(nb_video)

	if nb_video > nb_video_in_folder:
		print("Error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	video_ids = random.sample(range(1, nb_video_in_folder), nb_video)

	for id in tqdm(video_ids):
		video = video_in_folder[id]
		image_folder_name = video.split('.')[0]
		image_folder_name = image_folder_name.split('/')[-1]
		image_folder_path = join(sub_path, image_folder_name)
		extract_I_frames(video, image_folder_path)



if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	p.add_argument('--nb_video', '-n', type=str, default='10')
	args = p.parse_args()

	process(**vars(args))