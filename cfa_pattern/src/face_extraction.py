#!/usr/bin/env python

'''
face detection using haar cascades

USAGE:
    facedetect.py [--cascade <cascade_fn>] [--nested-cascade <cascade_fn>] [<video_source>]
'''

# Python 2/3 compatibility
from __future__ import print_function

import numpy as np
import cv2 as cv

# local modules
from video import create_capture
from common import clock, draw_str


def detect(img, cascade):
    rects = cascade.detectMultiScale(img, scaleFactor=1.3, minNeighbors=4, minSize=(30, 30),
                                     flags=cv.CASCADE_SCALE_IMAGE)
    if len(rects) == 0:
        return []
    rects[:,2:] += rects[:,:2]
    return rects

#def draw_rects(img, rects, color):
#    for x1, y1, x2, y2 in rects:
#        cv.rectangle(img, (x1, y1), (x2, y2), color, 2)

def face_extraction(args, video_src, faces):
    import sys, getopt

    try:
        video_src = video_src[0]
    except:
        video_src = 0
    args = dict(args)
    cascade_fn = args.get('--cascade', "data/haarcascades/haarcascade_frontalface_alt.xml")
    nested_fn  = args.get('--nested-cascade', "data/haarcascades/haarcascade_eye.xml")

    cascade = cv.CascadeClassifier(cv.samples.findFile(cascade_fn))
    nested = cv.CascadeClassifier(cv.samples.findFile(nested_fn))

    cam = create_capture(video_src, fallback='synth:bg={}:noise=0.05'.format(cv.samples.findFile('data/lena.jpg')))
    _ret, img = cam.read()
    while _ret:

        #convert to gray scale
        gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
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
                if len(subrects) == 0:
                    #store potential face as detected face
                    faces.put(img[y1:y2, x1:x2])


        _ret, img = cam.read()

    print('Faces detection - done')

if __name__ == '__main__':
    print("face_extraction stand alone")
    import sys, getopt
    args, video_src = getopt.getopt(sys.argv[1:], '', ['cascade=', 'nested-cascade='])
    #face_extraction(args, video_src)
    cv.destroyAllWindows()