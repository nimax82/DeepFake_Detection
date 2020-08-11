#LANGLADE Maxime
#26/02/20
import os
import numpy as np
import dlib
import cv2
from os.path import join
import argparse
from tqdm import tqdm
import glob

def count_image(list_folders):
	generic_image_name = "*.png"
	mean = 0
	for folder in tqdm(list_folders):
		#images_list = [f for f in glob.glob(join(folder, generic_image_name)) if "0*" == f or "1*" == f or "2*" == f or "3*" == f or "4*" == f or "5*" == f or "6*" == f or "7*" == f or "8*" == f]
		images_list = glob.glob(join(folder, generic_image_name))
		mean = mean + len(images_list)

	return mean / len(list_folders)


def process():
	p_folder = "/Volumes/VERBATIM_HD/Stage_Maxime/faceForensics/dataSet_c23/gen_test/images/P"
	i_folder = "/Volumes/VERBATIM_HD/Stage_Maxime/faceForensics/dataSet_c23/gen_test/images/I"
	b_folder = "/Volumes/VERBATIM_HD/Stage_Maxime/faceForensics/dataSet_c23/gen_test/images/B"

	generic_video_folder = "gen_*"
	
	list_folders_p = glob.glob(join(p_folder, generic_video_folder))
	list_folders_i = glob.glob(join(i_folder, generic_video_folder))
	list_folders_b = glob.glob(join(b_folder, generic_video_folder))
	
	print("count p:")
	list_cout_image_p = count_image(list_folders_p)
	print(list_cout_image_p)
	print("count i:")
	list_cout_image_i = count_image(list_folders_i)
	print(list_cout_image_i)
	print("count b:")
	list_cout_image_b = count_image(list_folders_b)
	print(list_cout_image_b)


if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	#p.add_argument('--data_path','-p' , type=str)
	args = p.parse_args()

	process(**vars(args))
