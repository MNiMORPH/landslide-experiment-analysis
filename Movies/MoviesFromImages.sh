#! /bin/sh

# Run inside sub-folder containing the images
ffmpeg -framerate 24 -i ImgSec%*.tif -c:v libx264 -vf scale=320:220 -pix_fmt yuv420p ../movie_images_24fps.mp4
