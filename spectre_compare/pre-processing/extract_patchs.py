#LANGLADE Maxime
#11/03/20
import os
from math import floor
from os.path import join
import argparse
import subprocess
import cv2
import random
import glob

from tqdm import tqdm

#global variables
frame_per_path = 20

#Crop frame frame in center in order to fit mask size
def create_patch(image, mask_w, mask_h):
	height, width, c = image.shape
	x_center = floor(1/2 * width)
	y_center = floor(1/2 * height)
	x_start = x_center - floor(1/2 * mask_w)
	y_start = y_center - floor(1/2 * mask_h)
	y_end = int(y_start + mask_h)
	x_end = int(x_start + mask_w)
	return image[int(y_start):y_end, int(x_start):x_end]

#Process frame extraction for the real video and this associate_fake
def extract_frames_from_couples(output_path, video_path, associate_fake, mask_w, mask_h):
	#print(data_path, " ", x_center, " " ,y_center)
	video_name = video_path.split("/")[-1]
	video_name = video_name.split(".")[0]
	#print(video_name)
	sub_output = join(output_path, video_name)
	os.makedirs(sub_output, exist_ok=True)

	os.makedirs(join(sub_output, 'real'), exist_ok=True)
	reader = cv2.VideoCapture(video_path)
	for i in range(frame_per_path):
		success, image = reader.read()		
		patch = create_patch(image, mask_w, mask_h)

		cv2.imwrite(join(sub_output, 'real', '{:04d}.tiff'.format(i)), patch)

	reader.release()

	os.makedirs(join(sub_output, 'synthesis'), exist_ok=True)
	for id_fake in associate_fake:
		video_fake_name = id_fake.split("/")[-1]
		video_fake_name = video_fake_name.split(".")[0]
		os.makedirs(join(sub_output, 'synthesis', video_fake_name), exist_ok=True)
		reader = cv2.VideoCapture(id_fake)
		for i in range(frame_per_path):
			success, image = reader.read()		
			patch = create_patch(image, mask_w, mask_h)
			height, width, c = patch.shape
			cv2.imwrite(join(sub_output, 'synthesis', video_fake_name, '{:04d}.tiff'.format(i)), patch)

		reader.release()

#search the smallest video resolution in order to crop all videos to this size
def get_smallest_resolution(folder_path_real, video_ids, fake_videos_paths):
	min_w = 1000
	min_h = 1000

	print("search for the smallest resolution")
	
	for i in video_ids:
		video = join(folder_path_real, os.listdir(folder_path_real)[i])
		reader = cv2.VideoCapture(video)
		success, image = reader.read()
		reader.release()
		if success:
			height, width, c = image.shape 
			min_w = min(min_w, width)
			min_h = min(min_h, height)

	for sub in fake_videos_paths:
		for j in sub:
			video = j
			reader = cv2.VideoCapture(j)
			success, image = reader.read()
			reader.release()
			if success:
				height, width, c = image.shape 
				min_w = min(min_w, width)
				min_h = min(min_h, height)

	print("Done (" , min_w, "-", min_h, ")")

	return min_w, min_h

'''
	nb_video = len(video_ids) + (len(fake_videos_paths) * len(fake_videos_paths[0]))

	for i in tqdm(range(nb_video)):
		if i < len(video_ids):
			video = join(folder_path_real, os.listdir(folder_path_real)[video_ids[i]])
		else:
			a = (i - len(video_ids)) % len(fake_videos_paths[0])
			b = (i - len(video_ids)) % len(fake_videos_paths)
			print(a , " - " , b)
			video = fake_videos_paths[a][b]
'''
		#print(video)

		

	#print(video_ids)
	#print(fake_videos_path)
	
#search associate (synthetise) videos from an original one
def get_couples(video_ids, videos_path_real, videos_path_synth, nb_compare):
	print("get couple videos:")
	associate_videos_container = [] #store associate videos for all real videos

	video_to_remove_ids = []

	for id_real in tqdm(video_ids):
		video_container = [] #store associate videos for current real video
		video = os.listdir(videos_path_real)[id_real]
		current_real_path = join(videos_path_real, video)

		name_video_fake = video.split('_')[0] + "_id*_" + video.split('_')[-1]
		#print(name_video_fake)
		potential_fake_videos = glob.glob(join(videos_path_synth, name_video_fake))

		if len(potential_fake_videos) < nb_compare:
			#print("error: not enough fake video for this real one")
			#video_ids.remove(id_real)
			video_to_remove_ids.append(id_real)
			continue

		#print(potential_fake_videos)

		fake_ids = random.sample(range(0, len(potential_fake_videos)), nb_compare)

		for i in fake_ids:
			video_container.append(potential_fake_videos[i])
		

		associate_videos_container.append(video_container)

	for id in video_to_remove_ids:
		video_ids.remove(id)

	return associate_videos_container	

	
#remove of the video list the ones who are not valid (unreadable)
def check_videos(video_ids, videos_path):
	print("select videos in folders:")

	for id in tqdm(range(len(video_ids))):

		off_set = 0
		success = False

		while not success:
			#open video (with offset if needed)
			id_r = (video_ids[id] + off_set) % len(os.listdir(videos_path))

			video = os.listdir(videos_path)[id_r] #todo check unicity in video_ids array

			off_set += 1

			# avoid hidden files
			if video[0] == '.':
				continue

			# try to read first frames
			reader = cv2.VideoCapture(join(videos_path, video))
			success, image = reader.read()
			reader.release()

			if not success:
				continue

			#reacheble only if video is correct
			video_ids[id] = id_r

		#todo check unicity in video_ids array

#main
def extract_videos(data_path, nb_video, nb_compare):
	folder_path_real = join(data_path, 'Celeb-real/videos')
	folder_path_synth = join(data_path, 'Celeb-synthesis/videos')
	images_path = join(data_path, 'image_couples')
	nb_compare = int(nb_compare)

	os.makedirs(images_path, exist_ok=True)

	nb_video_in_folder = len(os.listdir(folder_path_real))
	nb_video = int(nb_video)

	if nb_video > nb_video_in_folder:
		print("Error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	video_ids = random.sample(range(1, nb_video_in_folder), nb_video)

	check_videos(video_ids, folder_path_real)

	print(len(video_ids))

	fake_videos_path = get_couples(video_ids, folder_path_real, folder_path_synth, nb_compare)

	mask_w, mask_h = get_smallest_resolution(folder_path_real, video_ids, fake_videos_path)

	print("frame extraction:")
	for id in tqdm(range(len(video_ids))):
		video_path = join(folder_path_real, os.listdir(folder_path_real)[video_ids[id]])
		extract_frames_from_couples(images_path, video_path, fake_videos_path[id], mask_w, mask_h)



if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--data_path','-p' , type=str)
	p.add_argument('--nb_video', '-n', type=str, default='30')
	p.add_argument('--nb_compare', '-c', type=str, default='1')
	args = p.parse_args()

	extract_videos(**vars(args))


