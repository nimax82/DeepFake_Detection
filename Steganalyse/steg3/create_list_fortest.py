import random
import glob
from os.path import join

def get_associate(path_orig):
	path_fake = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/videos'
	name_orig = path_orig.split('/')[-1]
	name_orig = name_orig.split('.')[0]

	#print(name_orig)

	generic_fake_name = name_orig.split('_')[0] + "_id*_*"

	fakes_candidates = glob.glob(join(path_fake, generic_fake_name))

	
	ids = random.sample(range(1, len(fakes_candidates)), 1)
	return fakes_candidates[ids[0]]


def process(nb_video):

	video_list_path = []

	path_orig = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/videos'
	#path_fake = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/videos'

	generic_video_name = "id*"
	videos_folders_orig = glob.glob(join(path_orig, generic_video_name))
	#videos_folders_fake = glob.glob(join(path_fake, generic_video_name))
	#nb_video_in_folder = len(videos_folders)

	list_video_training_orig_init = glob.glob(join('/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/training', "id*"))
	list_video_training_orig = []
	for name in list_video_training_orig_init:
		if name.count("id") == 1:
			name_format = name.split('/')[-1]
			list_video_training_orig.append(name_format)

	orig_videos = []
	for potential_video in videos_folders_orig:
		name_video = potential_video.split('/')[-1]
		name_video = name_video.split('.')[0]
		if name_video not in list_video_training_orig:
			orig_videos.append(potential_video)

	#print(orig_videos)

	video_ids_orig = random.sample(range(1, len(orig_videos)), int(nb_video / 2))
	



	for id in video_ids_orig:

		video_list_path.append(orig_videos[id])
		video_list_path.append(get_associate(orig_videos[id]))


		#video_list_path.append(fake_videos[videos_ids_fake[i]])


	return video_list_path






