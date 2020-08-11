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
	select_cmd = " -f image2 -vf " + '''"select='eq(pict_type,PICT_TYPE_P)'"'''
	#select_cmd = " -f image2 -vf " + '''"select='eq(pict_type,B)'"'''
	cmd = "ffmpeg -i " + video_path + select_cmd + " -vsync vfr " + out_file
	
	os.system(cmd)              


def process(data_path):
	#video_list_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/test'
	#video_list = glob.glob(join(video_list_path, 'id*'))

	#video_to_process = []
	#for video_path in video_list:
	#	name = video_path.split('.')[0]
	#	name = name.split('/')[-1]
	#	video_to_process.append(name)

	video_path = join(data_path, "videos")
	images_path = join(data_path, "images", "resize_unl_P")

	video_list = glob.glob(join(video_path, 'gen_*'))
	#video_list = [video_list_orig, video_list_fake]

	print(video_list)


	os.makedirs(images_path, exist_ok=True)


	#video_in_folder = glob.glob(join(videos_path, generic_video_name))

	#nb_video_in_folder = len(video_to_process)

	for video in tqdm(video_list):
		video_name = video.split('.')[0]
		video_name = video_name.split('/')[-1]
		output_folder = join(images_path, video_name)
		extract_I_frames(video, output_folder)


if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	args = p.parse_args()

	process(**vars(args))