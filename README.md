# landslide-experiment-analysis

Codes written to manage, analyze, and display data from autogenic landslides in the 2016--2017 sandbox experiments at SAFL.

* **Georeference**: Code to batch-georeference the orthophotos to the DEMs based on a set of hand-picked control points. To gather these, we used the QGIS georeferencing tool.

* **GrainSize**: Contains and plots data on the grain-size distribution used in the experiment.

* **Landslides**: Code to analyze and plot landslide geometries, chronologies, and derived relationships. This folder contains its own README with information on its contents.

* **Movies**: *ffmpeg* commands to build movies of the experiment from the images; to be run in the folders containing the images, and will save the output in the next directory up.

**Chronology**: For the key scripts needed used to semi-automate the generation of continuous and consistent time-series data from disparate experiment measurements and measurement types, see https://github.com/umn-earth-surface/fluvial-experiment-time-series-helper. This latter version is the version of record.

