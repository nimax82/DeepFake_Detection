# Python 2/3 compatibility
from __future__ import print_function

import numpy as np
import cv2 as cv

import sys, getopt

# local modules
from thrid.video import create_capture
from thrid.common import clock, draw_str


def detect(img, cascade):
	rects = cascade.detectMultiScale(img, scaleFactor=1.3, minNeighbors=4, minSize=(30, 30),
									 flags=cv.CASCADE_SCALE_IMAGE)
	if len(rects) == 0:
		return []
	rects[:,2:] += rects[:,:2]
	return rects

def draw_rects(img, rects, color):
	for x1, y1, x2, y2 in rects:
		cv.rectangle(img, (x1, y1), (x2, y2), color, 2)

def face_extraction(frame):
	

	cascade_fn = "data/haarcascades/haarcascade_frontalface_alt.xml"
	nested_fn  = "data/haarcascades/haarcascade_eye.xml"

	cascade = cv.CascadeClassifier(cv.samples.findFile(cascade_fn))
	nested = cv.CascadeClassifier(cv.samples.findFile(nested_fn))

	#cam = create_capture(video_src, fallback='synth:bg={}:noise=0.05'.format(cv.samples.findFile('pre-processing/data/lena.jpg')))
	#_ret, img = cam.read()
	



	#convert to gray scale
	gray = cv.cvtColor(frame, cv.COLOR_BGR2GRAY)
	gray = cv.equalizeHist(gray)

	# detected faces
	rects = detect(gray, cascade)

	#In order to avoid false detection, need to check if the potential faces have eyes inside.
	vis = img.copy()
	if not nested.empty():
		for x1, y1, x2, y2 in rects:
			roi = gray[y1:y2, x1:x2]
			vis_roi = vis[y1:y2, x1:x2]
			subrects = detect(roi.copy(), nested)
			if len(subrects) >= 1:

			#x_center = int((x2 + x1) / 2)
			#y_center = int((y2 + y1) / 2)
				return x1, y1, x2, y2

	return -1,-1
	print('Error: no face detected')

if __name__ == '__main__':
	cv.destroyAllWindows()