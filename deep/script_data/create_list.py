import random
import glob
import cv2
from os.path import join

def get_associate(path_orig, is_multi):
	path_fake = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/videos'
	name_orig = path_orig.split('/')[-1]
	name_orig = name_orig.split('.')[0]

	#print(name_orig)

	generic_fake_name = name_orig.split('_')[0] + "_id*_*"

	fakes_candidates = glob.glob(join(path_fake, generic_fake_name))

	#fake_to_return = []

	print(len(fakes_candidates))

	#ids = random.sample(range(0, len(fakes_candidates)), len(fakes_candidates))

	count_id = 0
	
	is_valid = False
	while not is_valid:
		reader = cv2.VideoCapture(fakes_candidates[count_id])
		success, image = reader.read()
		w, h, c = image.shape
		if w >= 299 and h >= 299:
			is_valid = True
		else:
			count_id += 1

	if not is_multi:
		return fakes_candidates[count_id]

	else:
		count_id_2 = count_id + 1
	
		is_valid = False
		while not is_valid:
			reader = cv2.VideoCapture(fakes_candidates[count_id_2])
			success, image = reader.read()
			w, h, c = image.shape
			if w >= 299 and h >= 299:
				is_valid = True
			else:
				count_id_2 += 1

		return fakes_candidates[count_id], fakes_candidates[count_id_2]


def process(nb_orig):

	video_list_path = []

	path_orig = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/videos'
	generic_video_name = "id*"
	videos_folders = glob.glob(join(path_orig, generic_video_name))
	nb_video_in_folder = len(videos_folders)

	video_ids = random.sample(range(1, nb_video_in_folder), nb_orig)

	for i in range(nb_orig):
		video_list_path.append(videos_folders[video_ids[i]])
		fake = get_associate(videos_folders[video_ids[i]], bool(i%2))
		if bool(i%2):
			video_list_path.append(fake[0])
			video_list_path.append(fake[1])
		else:
			video_list_path.append(fake)

	return video_list_path






