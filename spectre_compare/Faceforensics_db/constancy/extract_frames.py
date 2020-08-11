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


def extract_frame(video):

	out_path = "/Volumes/VERBATIM HD/Stage_Maxime/faceForensics/picked_for_fft/constancy/images"

	
	reader = cv2.VideoCapture(video)
	frame_num = 0

	while reader.isOpened():
		success, image = reader.read()
		frame_num += 1
		if not success:
			#print("frame error for: ", data_path, ' frame n° ', frame_num)
			break

		cv2.imwrite(join(out_path, '{:04d}.tiff'.format(frame_num)), image)

	
	

	#print(id_video)
	#export_to_json(["eye_distance","area"], [eye_distances, face_areas], id_video, id_mode)
	#export_to_json("eye_distance", eye_distances)
	#export_to_json("area", face_areas)		



def extract_frames(video_path):

	
	extract_frame(video_path)




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--video_path','-p' , type=str)
	args = p.parse_args()

	extract_frames(**vars(args))


