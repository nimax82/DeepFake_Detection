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


import json


selected_frame = 70
padding = 50

def rect_to_bb(rect):
	# take a bounding predicted by dlib and convert it
	# to the format (x, y, w, h) as we would normally do
	# with OpenCV
	x = rect.left()
	y = rect.top()
	w = rect.right() - x
	h = rect.bottom() - y
	# return a tuple of (x, y, w, h)
	return (x, y, w, h)


def shape_to_np(shape, dtype="int"):
	# initialize the list of (x, y)-coordinates
	coords = np.zeros((68, 2), dtype=dtype)
	# loop over the 68 facial landmarks and convert them
	# to a 2-tuple of (x, y)-coordinates
	for i in range(0, 68):
		coords[i] = (shape.part(i).x, shape.part(i).y)
	# return the list of (x, y)-coordinates
	return coords


def crop_and_save_frame(output_path, image, bbox, out_name):

	x_start = bbox[0] - padding
	x_end = x_start + bbox[2] + (2 * padding)
	y_start = bbox[1] - padding
	y_end = y_start + bbox[3] + (2 * padding)
	patch = image[y_start:y_end, x_start:x_end]
	cv2.imwrite(join(output_path, out_name + '.tiff'), patch)





def extract_frame(orig_path, video_orig, fake_path, video_fake):

	# initialize dlib's face detector (HOG-based) and then create
	# the facial landmark predictor
	detector = dlib.get_frontal_face_detector()
	predictor = dlib.shape_predictor("data/shape_predictor_68_face_landmarks.dat")

	reader_orig = cv2.VideoCapture(video_orig)
	reader_fake = cv2.VideoCapture(video_fake)
	frame_num = 0

	rects_orig = []
	rects_fake = []

	while (reader_orig.isOpened() and reader_fake.isOpened() and frame_num < selected_frame - 1):
		success_o, image_fake = reader_orig.read()
		success_f, image_orig = reader_fake.read()
		frame_num += 1
		success = success_o | success_f

		if not success:
			#print("frame error for: ", data_path, ' frame n° ', frame_num)
			break

	


	while len(rects_orig) == 0 and len(rects_fake) == 0:
		success_o, image_orig = reader_orig.read()
		success_f, image_fake = reader_fake.read()

		gray_orig = cv2.cvtColor(image_orig, cv2.COLOR_BGR2GRAY)
		gray_fake = cv2.cvtColor(image_fake, cv2.COLOR_BGR2GRAY)

		# detect faces in the grayscale image
		rects_orig = detector(gray_orig, 1)
		rects_fake = detector(gray_fake, 1)

	reader_orig.release()
	reader_fake.release()

	bbox_orig = rect_to_bb(rects_orig[0])
	bbox_fake = rect_to_bb(rects_fake[0])

	out_name = video_orig.split("/")[-1][0]
	crop_and_save_frame(orig_path, image_orig, bbox_orig, out_name)
	out_name = video_fake.split("/")[-1][0]
	crop_and_save_frame(fake_path, image_fake, bbox_fake, out_name)
	

	#print(id_video)
	#export_to_json(["eye_distance","area"], [eye_distances, face_areas], id_video, id_mode)
	#export_to_json("eye_distance", eye_distances)
	#export_to_json("area", face_areas)		


def get_associate(input_video_path, fake_folder_path):
	input_video_name = input_video_path.split('/')[-1]
	id_input = input_video_name.split('_')[0] + "_*"
	fake_video_name = glob.glob(join(fake_folder_path, id_input))
	return fake_video_name[0]


def extract_frames(videos_path):

	orig_path = join(videos_path, "orig")
	fake_path = join(videos_path, "fake")

	videos = [vid for vid in glob.glob(join(orig_path, "*_0*")) if ".*.mp4" not in vid]
	
	for video_orig in tqdm(videos):
		video_fake = get_associate(video_orig, fake_path)
		extract_frame(orig_path, video_orig, fake_path, video_fake)




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--videos_path','-p' , type=str)
	args = p.parse_args()

	extract_frames(**vars(args))


