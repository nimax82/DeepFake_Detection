#LANGLADE Maxime
#17/06/20
import os
from os.path import join
from shutil import copyfile
import argparse
import glob

nb_I = 3
nb_P = 25
input_folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/faceForensics/dataSet_c23/gen_test/images'
output_folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/faceForensics/dataSet_c23/gen_test/images/resize'

def copy_frames(folder_path):
	folder_name = folder_path.split("/")[-1]
	output_path = join(output_folder_path, folder_name)
	os.makedirs(output_path, exist_ok=True)

	generic_image_name = "0*.png"
	images_list = glob.glob(join(folder_path, generic_image_name))

	counter_image = 0

	for image in images_list:
		if counter_image == nb_I:
			break
		#print("current image path" , image)

		image_name = image.split("/")[-1]
		image_name = 'fi_' + image_name

		#print("copy: ", image_name)
		image_dst = join(output_path, image_name)
		#print("to: ", image_dst)

		copyfile(image, image_dst)

		counter_image += 1
 	
	folder_P = join(input_folder_path, 'P', folder_name)
	images_list = glob.glob(join(folder_P, generic_image_name))

	counter_image = 0

	for image in images_list:
		if counter_image == nb_P:
			break
		#print("current image path" , image)

		image_name = image.split("/")[-1]
		image_name = 'fp_' + image_name
		#print("copy: ", image_name)
		image_dst = join(output_path, image_name)
		#print("to: ", image_dst)

		copyfile(image, image_dst)

		counter_image += 1



def process():
	generic_video_name = "gen_*"
	images_folders = glob.glob(join(input_folder_path, 'I',generic_video_name))

	for id in images_folders:
		copy_frames(id)




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	#p.add_argument('--data_path','-p' , type=str)
	#p.add_argument('--nb_video', '-n', type=str, default='50')
	args = p.parse_args()

	process(**vars(args))
