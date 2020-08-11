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
import create_list_fortest

from tqdm import tqdm            


def process(data_path, nb_video):
	nb_video = int(nb_video)
	video_to_process = create_list_fortest.process(nb_video)
	


	#sub_path = join(data_path, "training")
	deep_folder = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/test'
	os.makedirs(deep_folder, exist_ok=True)


	
	#video_in_folder = glob.glob(join(videos_path, generic_video_name))

	#nb_video_in_folder = len(video_to_process)

	for video in tqdm(video_to_process):
		image_folder_name = video.split('.')[0]
		image_folder_name = image_folder_name.split('/')[-1]
		image_folder_path = join(sub_path, image_folder_name)



if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	p.add_argument('--nb_video', '-n', type=str, default='100')
	args = p.parse_args()

	process(**vars(args))