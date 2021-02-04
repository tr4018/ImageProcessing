# A1-Image-Processing
Talia Rahall, Tom Andrews A1 Image Processing Jan 2021

Code dependencies: Astopy.io (photoutils) and pandas

A1_Main.py
Contains the main program, this can simply be run and will scan image, identify objects and perform analysis.

A1_global_background.py
Produces a histogram for the large pixel counts of the image
Produces a histogram of the background noise and calculates the mean value using a gaussian fit
Recomended that you walk through the code in cells

A1_testing.py
Tests the object_find dunction a 2d synthetic gaussian (galaxy), and determines effective radius of synthetic galaxy
Displays testing on a small region of the image, displaying our aperture photometry method operating successfully on two adjacent galaxies
Recomended that you walk through the code in cells, particularly when inputting a your own box to scan, therefore you can see the aperture photometry operate on each identified object with the produced plots.
