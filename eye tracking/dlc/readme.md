## a few scripts to help with using [Deep Lab Cut](https://github.com/DeepLabCut/DeepLabCut) for eye tracking

- improveContrast.m : boost the contrast of eye movies before DLC to facilitate manual labeling step 
- rdk_2afc_eye_tracking.ipynb : based on one of the example DLC files for training a model and labeling new videos. 
  fits 8 points to the border of the pupil (need >=5 to fit ellipse), 2 to the corneal reflection and 1 in each eye corner
  ![like so](https://github.com/lpetreanu/petreanulab/blob/master/eye%20tracking/dlc/example%20labeled%20images/Training-ec_MF122_180215_eye_4b1-img768.png?raw=true)
- config.yaml : an example config file for the model
- read_dlc.m : matlab script to read the results of DLC (saved as .csv) and fit an ellipse
