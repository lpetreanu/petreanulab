### Here's how to use this set of scripts to transform your video recordings of the mouse's eye into eye tracking data!

# The algorithm for eye tracking used in this code comes from Zoccolan et al Frontiers paper.
# Check out this figure to have an understanding of how it works: 
# https://www.frontiersin.org/files/Articles/2035/fnins-04-00193-HTML/image_m/fnins-04-00193-g002.jpg



# Step 1: crop your eye videos to only include the animal's eye
# eyeVidChopper.m does this and aligns to frames from ScanImage (if you use SI to trigger the eye camera frames)
# eyePrep.m does this for other eye movies. you can input all the sessions you have and it will crop them all 
# (script is currently from Gabi's FB silencing experiment, but can be adapted for other experiments)

# Step 2: figure 
