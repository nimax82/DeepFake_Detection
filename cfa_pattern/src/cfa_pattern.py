import threading
import cv2 as cv
import queue

from face_extraction import face_extraction	

#Extract all faces from the input video and return them as list of cv:Mat
def buildFaces():
	import sys, getopt
	args, video_src = getopt.getopt(sys.argv[1:], '', ['cascade=', 'nested-cascade='])
	faces_queue = queue.Queue()
	face_extraction_process = threading.Thread(target=face_extraction, args=(args, video_src, faces_queue))
	face_extraction_process.start()
	print("started")
	
	faces = []

	#todo delete join...
	face_extraction_process.join()
	isAllAdded = False
	while not isAllAdded:
		print(faces_queue.qsize())
		faces.append(faces_queue.get())
		if faces_queue.qsize() > 0:
			isAllAdded = False 
		elif not face_extraction_process.isAlive():
			isAllAdded = True
		
		faces_queue.task_done()
		#cv.imshow('detected faces', faces[-1]) #last
		#print("press any key...")
		#cv.waitKey(0)

	#face_extraction_process.join()
	print("all faces added")
	return faces

def main():
	faces = buildFaces()
	#cfa...


if __name__ == '__main__':
    main()
