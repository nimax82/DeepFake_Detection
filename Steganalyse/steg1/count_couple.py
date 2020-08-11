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


def count_type(video):
	global count_alone
	global count_couple
	#video_name = video.split(".")[0]
	associate_real = video_name.split("_")[0] + "_" + video_name.split("_")[-1]
	print(associate_real)
	path_training = join("/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images/training", associate_real)
	print(path_training)
	#training = glob.glob(join("/Volumes/VERBATIM\ HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images/training", associate_real))
	training = os.path.isdir(path_training)
	print(training)
	if training:
		if count_couple >= 50:
			return False
		print("couple")
		count_couple += 1
	else:
		if count_alone >= 50:
			return False
		print("alone")
		count_alone += 1

	return True


def extract_videos(data_path):
	nb_video_in_folder = len(os.listdir(data_path))

	for id in tqdm(range(len(video_ids))):
		video_name = video_ids[id].split("/")[-1]
		if video_name == '.':
			continue
		print(video_name)
		exit()
		count_type(video_name)
		

	print("couple = ", count_couple)
	print("alone = ", count_alone)



if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	args = p.parse_args()

	extract_videos(**vars(args))


