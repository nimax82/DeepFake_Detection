#LANGLADE Maxime
#26/02/20
import os
from os.path import join
from shutil import copyfile
import argparse
from tqdm import tqdm
import glob

nb_frame = 15
input_folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/test'
output_folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/limited_test'

def copy_frames(folder_path):
	folder_name = folder_path.split("/")[-1]
	output_path = join(output_folder_path, folder_name)
	os.makedirs(output_path, exist_ok=True)

	generic_image_name = "0*.tiff"
	images_list = glob.glob(join(folder_path, generic_image_name))

	counter_image = 0

	for image in images_list:
		if counter_image == nb_frame:
			return
		#print("current image path" , image)

		image_name = image.split("/")[-1]
		#print("copy: ", image_name)
		image_dst = join(output_path, image_name)
		#print("to: ", image_dst)

		copyfile(image, image_dst)

		counter_image += 1



def process():
	generic_video_name = "id*"
	images_folders = glob.glob(join(input_folder_path, generic_video_name))

	for id in tqdm(images_folders):
		#print("new folder image")
		copy_frames(id)




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	#p.add_argument('--data_path','-p' , type=str)
	#p.add_argument('--nb_video', '-n', type=str, default='50')
	args = p.parse_args()

	process(**vars(args))
