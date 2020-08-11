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
	video_list_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/test'
	video_list = glob.glob(join(video_list_path, 'id*'))

	video_to_process = []
	for video_path in video_list:
		name = video_path.split('.')[0]
		name = name.split('/')[-1]
		video_to_process.append(name)


	sub_path = join('/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/I_P_Compare/P', "test")
	#sub_path = join(images_path, 'testV2')
	os.makedirs(sub_path, exist_ok=True)


	#video_in_folder = glob.glob(join(videos_path, generic_video_name))

	#nb_video_in_folder = len(video_to_process)

	for video in tqdm(video_to_process):
		if video.count("id") == 2:
			sub_folder = "Celeb-synthesis/videos"
		else:
			sub_folder = "Celeb-real/videos"
		image_folder_path = join(sub_path, video)
		video_ext = video + ".mp4"
		extract_I_frames(join(data_path, sub_folder , video_ext), image_folder_path)



if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	args = p.parse_args()

	process(**vars(args))