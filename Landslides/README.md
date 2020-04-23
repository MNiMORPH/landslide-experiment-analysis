* `landslide_geometry.py` should be run in subfolders containing each of the shapefile sets for each experiment. It generates comma-separated files containing information on the position, shape (using `ellipses.py` among other approaches), orientation, and timing of each landslide event.
* `ellipses.py` fits an ellipse to each set of (x, y) points.
* `all_landslides_stats.py` should be run next, and from the folder one level above the indivudual-experiment-run subfolders (each of which will store a `Landslides.txt` output file with the comma-separated outputs produced above) to generate sets of statistics and plots, including those used in the main text.
* `width_area.py` relates landslide width and area; it an extra plotting script that is very similar to `all_landslides_stats.py`.
* `Nlandslides.py` simply counts the numbers of landslides in each experiment.
