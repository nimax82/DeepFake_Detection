#LANGLADE Maxime
#26/02/20
import os
import numpy as np
import dlib
import cv2
from os.path import join
import argparse
import random

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

def get_ldmarks():

	# initialize dlib's face detector (HOG-based) and then create
	# the facial landmark predictor
	detector = dlib.get_frontal_face_detector()
	predictor = dlib.shape_predictor("data/shape_predictor_68_face_landmarks.dat")

	reader = cv2.VideoCapture(0)
	while (reader.isOpened()):
		success, image = reader.read()

		#image = imutils.resize(image, width=500)
		gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

		# detect faces in the grayscale image
		rects = detector(gray, 1)
		#x1, y1, x2, y2 = face_extraction(frame):

		if len(rects) > 0:
			shape = predictor(gray, rects[0])
			shape = shape_to_np(shape)

			(x, y, w, h) = rect_to_bb(rects[0])

			# loop over the (x, y)-coordinates for the facial landmarks
			# and draw them on the image
			for (x, y) in shape:
				cv2.circle(image, (x, y), 1, (0, 0, 255), -1)
			

			#show the output image with the face detections + facial landmarks
			cv2.imshow("Output", image)
			cv2.waitKey(1)
	
	
	reader.release()		



if __name__ == '__main__':
	get_ldmarks()


