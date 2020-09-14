Analysis pipeline (as of June 15, 2020)

1. Run Suite2p
	ops file in C:\Users\Camille Mazo\suite2p
	classifier in E:\TempData
2. Run "AfterSuite2p_v2.m". Replaces "AfterSuite2p_singlePlane.m"
	/!\ Don't delete FlyBack frame folder, because the code is discarding the last plane
	Works for both single plane and multiple planes. If multiple planes, assumes that their is a FlyBack frame to be discarded.
	Input: Solely change the main folder name and "saveName"
	Output: Ftraces_all
	Also generates a log file
3. Run "SortData_PhotoStim_CM_v3.m". Modified for multiple planes (16 June 2020)
	Input: fTraces, SessionData, selection (1: only curated ROIs from s2p)
	Ouput: data. FileName is AnimalName_Position#
	Also plot a few Sanity Check figures
	And generates a log file

4a "NewAlignment_test.m"
	Will replace "zStack_avgPerPlane.m" if everything works well.
	/!\ removing bad plane worked with CMloop34_pos2  (pl.6) but not with CMloop34_pos4 (pl.2). threshold needs to be adapted to each image?
	a. Average frames per planes.
	b. Use measure of crispness to throw away bad planes.
	c. Compute max projection
	d. Register individual frame to the max projection
	Inputs: ischannel2 = 1 if there is a 2nd channel; O otherwise
	Outputs: path (to be used by "RegisterStackToFunctionalImaging_v2")
	(4. zStack_avgPerPlane deprecated
		Inputs: ischannel2 = 1 if there is a 2nd channel; O otherwise
			registration. Do FFT-based registration or not. Usually it's done in the next step.
		Outputs: path. It's a global and will be used directly as an input of the next function
			 plane-averaged z-stack, with channels splitted if ischannel2 = 1.)
--- OR -----
4b "fast_zStack_avg.m" if structural stack obtained using fast-z (to get the same tilt as during imaging)

5. RegisterStackToFunctionalImaging_v2
	/!\ Check the file names next time runing it after "NewAlignment_test.m"
	Inputs: data (output from "SortData_PhotoStim_CM_v3")
		crop. Number of pixels to crop from the left and right-end sides. Typically [0 0] for L2/3, [25 25] is good for L5 (to eliminate pixels with light contamination)
		regStack. Register all planes in the z-stack and save it / only save the best plane if set to 0.
		path. Path to the z-stacks
	Outputs: Tiff files. Registered best planes and average with one plane above/below + registred z-stack if regStack = 1.
		 saves the "rescale_factor" in data.s2p
	If the z-stack and the functional imaging are not acquired using the same pixel resolution (e.g, 512x512 for z-stack and 1024x1024 for imaging), it rescales automatically the functional imaging mean image.
	If the variable "path" is not in the workspace, the function asks where is the path to the stacks.
	Also Runs "PlotROIsMasks_2"
# 6. PlotROIsMasks_2

####### Manual Identification of bead pos cells ################

7. "IDfying_redCells.m" Also updated for multiple planes
	Inputs: data.mat file (output from 3.)
		x and y coordinates of the identified bead+ cells
	Outputs: saves data.beads_pos
		 saves the tiff of s2p meanImg overlaid with masks and cross over bead pos cells

8.2 "BestResp_BeadsSelec_v5.m"
	Inputs: data, selection, RespLag, comp_mthd, alpha, IndivROI, saveFig, SessionData
		selection is for the bead+ cells
		alpha is a vector with 1) threshold for vis resp (up to 4 paired t-tests) and 2) threshold for effect of opto on best significant vis resp (single paired t-test)
		IndivROI: Runs "PlotIndivROIs.m" with same options and IndivROIFigFlag =1, TuningCurveFigFlag =1
		saveFig for saving the quantifiaction plots
		SessionData only necessary if sfvalues ans tfvalues not save in data strucutre (version of SortData_PhotoStim_CM_v3 before April 15, 2020)
	Ouputs: 'quantif', 'h_vis' and 'dir_tuning', 'fit' if IndivROI =1. All saved in "data" structure
9. TuningData = TuningAnalysis_v2(TuningData,0,1)

####### Population Analysis ################
10. Run pop_quantif (for "ImpactOfLaser" figures)
	here is the threshold for response amplitude
	Inputs: data_list, saveFig
		data_list: cell array with the "quantif" (output of BestResp_BeadsSelec_v4 or later) from different sessions.
	
11. TuningData = pop_TuningAnalysis(data_list,do_MF,saveFig)
