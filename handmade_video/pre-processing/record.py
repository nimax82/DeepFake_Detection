import cv2

'''
vid_capture = cv2.VideoCapture(1)
vid_cod = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
output = cv2.VideoWriter("recording.mp4", vid_cod, 24.0, (1028,720))
'''
'''
fps = 15
capSize = (1028,720) # this is the size of my source video
fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v') # note the lower case
self.vout = cv2.VideoWriter()
success = self.vout.open('output.mov',fourcc,fps,capSize,True)
'''

# Create a video capture object and set some parameters
vid_capture = cv2.VideoCapture(1)
fps = 30.0
capsize = (1280, 720)

# Define the codec and create VideoWriter Object
fourcc = cv2.VideoWriter_fourcc('m', 'p','4','v')
output = cv2.VideoWriter()
success = output.open('output.mov', fourcc, fps, capsize, True)

while(success):
     # Capture each frame of webcam video
     success,frame = vid_capture.read()
     cv2.imshow("My cam video", frame)
     output.write(frame)
     # Close and break the loop after pressing "x" key
     if cv2.waitKey(1) &0XFF == ord('x'):
         break


# close the already opened camera
vid_capture.release()
# close the already opened file
output.release()
# close the window and de-allocate any associated memory usage
cv2.destroyAllWindows()