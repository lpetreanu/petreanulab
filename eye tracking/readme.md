## Here's how to use this set of scripts to transform your video recordings of the mouse's eye into eye tracking data!

The algorithm for eye tracking comes from 
> Zoccolan, D., Graham, B. J., & Cox, D. D. (2010). A self-calibrating, camera-based eye tracker for the recording of rodent eye movements. Frontiers in Neuroscience, 4(NOV), 1–12. https://doi.org/10.3389/fnins.2010.00193 

Check out [this figure](https://www.frontiersin.org/files/Articles/2035/fnins-04-00193-HTML/image_m/fnins-04-00193-g002.jpg) to understand how it works.

- EyeBatchAnalyzer.m is the file that will loop through all the videos you select and run the eye tracker on them. It will save a structure with eye tracking info for that session.
- analyzeThatPupil.m is the code that actually implements the eye tracking algorithm. It will save a video with the fitted pupil so you can see how well it did.

### Step 1: crop your eye videos to only include the animal's eye
- eyeVidChopper.m does this and aligns to frames from ScanImage (if you use SI to trigger the eye camera frames)
- eyePrep.m does this for other eye movies. you can input all the sessions you have and it will crop them all (script is currently from Gabi's FB silencing experiment, but can be adapted for other experiments). It will also fix the trial number.

### Step 2: adapt parameters to your lighting and recording conditions 

In EyeBatchAnalyzer, you will create a set of parameters for the eye tracking algorithm to use. This will be saved in the same folder as the eye movies, so you will always know what you used. Set the displayFlag=2 in EyeBatchAnalyzer, and the dbug =1 in analyzeThatPupil to visualize how the parameters affect your tracker.
_add explanations of parameters here_

### Step 3: visualize your results, detect saccades
This part is not very developed, but you can see some examples and first steps in eyePlot.m (agian, written for Gabi's FB silencing project, but can be adapted).


