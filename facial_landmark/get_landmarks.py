#LANGLADE Maxime
#26/02/20
import os
import numpy as np
import dlib
import cv2
from os.path import join
import argparse
import random
from tqdm import tqdm


from shapely.geometry import Polygon

import json

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

def export_to_json(key_list, data_list, id_video, mode):
	dict_to_export = {}
	for i in range(len(key_list)):
		dict_to_export[key_list[i]] = data_list[i]

	out_name = "./result/" + mode + "/" + id_video + ".json"
	with open(out_name, 'w') as json_file:
  		json.dump(dict_to_export, json_file)



#approx area of a list of points
def approx_area(shape):
	#bounding points are in range [0:26]
	edge_points = shape[:27]
	return Polygon(edge_points).area
	#print("area poly = ", poly.area)


#Caution: return only the X value
def convert_to_local(width, x_min, pos):
	(pos_X, pos_Y) = pos
	return (pos_X - x_min) / width


def compute_eyes_distances(shape):

	(edge_leftX, edge_leftY) = shape[16]
	(edge_rightX, edge_rightY) = shape[0]
	width_face =  edge_leftX - edge_rightX

	#36-39
	pupil_rr = convert_to_local(width_face ,edge_rightX, shape[36])
	pupil_rl = convert_to_local(width_face ,edge_rightX, shape[39])

	#42-45
	pupil_lr = convert_to_local(width_face ,edge_rightX, shape[42])
	pupil_ll = convert_to_local(width_face ,edge_rightX, shape[45])


	#left and right of the character
	pupil_left = (pupil_lr + pupil_ll) / 2
	pupil_right = (pupil_rr + pupil_rl) / 2


	return abs(pupil_right - pupil_left)

'''

	print("edge right global: ",  edge_rightX)
	print("edge right local: ",  convert_to_local(width_face ,edge_rightX, shape[0]))

	print("edge left global: ",  edge_leftX)
	print("edge left local: ",  convert_to_local(width_face ,edge_rightX, shape[16]))

	print("width_face: ", width_face)

	(p_rightRX, p_rightRY)  = shape[36]
	(p_rightLX, p_rightLY)  = shape[39]
	p_r_global = (p_rightLX + p_rightRX) / 2

	#print(p_rightRX)
	#print(p_rightLX)



	print("-----")
	print("p right global: " ,  p_r_global)
	print("p right local: " , pupil_right)

	print("p left: " , pupil_left)
'''




def get_ldmarks(data_path):

	# initialize dlib's face detector (HOG-based) and then create
	# the facial landmark predictor
	detector = dlib.get_frontal_face_detector()
	predictor = dlib.shape_predictor("data/shape_predictor_68_face_landmarks.dat")

	reader = cv2.VideoCapture(data_path)
	frame_num = 0

	eye_distances = []
	face_areas = []

	while (reader.isOpened()):
		success, image = reader.read()
		frame_num += 1

		if not success:
			#print("frame error for: ", data_path, ' frame n° ', frame_num)
			break

		#image = imutils.resize(image, width=500)
		gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

		# detect faces in the grayscale image
		rects = detector(gray, 1)
		#x1, y1, x2, y2 = face_extraction(frame):

		if len(rects) == 0:
			continue

		shape = predictor(gray, rects[0])
		shape = shape_to_np(shape)

		(x, y, w, h) = rect_to_bb(rects[0])

		# loop over the (x, y)-coordinates for the facial landmarks
		# and draw them on the image
		id_point = 0

		#print("eyes distance = ", compute_eyes_distances(shape))
		
		eye_distances.append(compute_eyes_distances(shape))
		face_areas.append(approx_area(shape))

		for (x, y) in shape:
			#print(x, "-" ,y)
			cv2.circle(image, (x, y), 1, (0, 0, 255), -1)
			#print(id_point)
			#id_point += 1
			#cv2.imshow("Output", image)
			#cv2.waitKey(0)
		

		#show the output image with the face detections + facial landmarks
		#cv2.imshow("Output", image)
		#cv2.waitKey(1)
			"""
		if cv2.waitKey(1) == ord('q'):
			reader.release()
			cv2.DestroyAllWindows()
			break
			"""
	
	
	reader.release()

	id_video = data_path.split('/')[-1]
	id_video = id_video.split('.')[0]

	id_mode = data_path.split('2/Celeb-')[-1]
	id_mode = id_mode.split('/videos')[0]

	print(id_mode)

	#print(id_video)
	export_to_json(["eye_distance","area"], [eye_distances, face_areas], id_video, id_mode)
	#export_to_json("eye_distance", eye_distances)
	#export_to_json("area", face_areas)		

def get_video(videos_path, id):
	off_set = 0
	success = False

	while not success:
		#print("video id: " + video)

		id_r = (id + off_set) % len(os.listdir(videos_path))
		video = os.listdir(videos_path)[id_r] #todo check unicity in video_ids array

		off_set += 1

		if video[0] == '.':
			continue

		reader = cv2.VideoCapture(join(videos_path, video))
		success, image = reader.read()
		reader.release()
		

	return video


def folder2landmarks(videos_path, nb_video):


	nb_video_in_folder = len(os.listdir(videos_path))
	nb_video = int(nb_video)

	if nb_video > nb_video_in_folder:
		print("Error: nb_video > nb_video_in_folder")
		nb_video = nb_video_in_folder

	video_ids = random.sample(range(1, nb_video_in_folder), nb_video)
	
	for id in tqdm(video_ids):
		video = get_video(videos_path , id)
		get_ldmarks(join(videos_path, video))




if __name__ == '__main__':
	p = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	p.add_argument('--videos_path','-p' , type=str)
	p.add_argument('--nb_video', '-n', type=str, default='1')
	args = p.parse_args()

	folder2landmarks(**vars(args))


