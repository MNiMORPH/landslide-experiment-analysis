#! /bin/sh

# Run inside sub-folder containing the images
ffmpeg -framerate 15 -i ImgSec%*.tif -c:v libx264 -vf scale=960:660 -pix_fmt yuv420p ../movie_images_960x660p_15fps__5minInExperiment_per_second.mp4
