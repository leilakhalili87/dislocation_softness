I.Adjust the X-coordinates to remove rigid-body translations (adjust_xcoor/)
	1. To remove the rigid-body displacement during the dislocation movement.
	2. Make directory in gb_id folder: mkdir adjusted_xcoor
	3. How to run: ovitos adjust_xcoor.py

II. Calculate PB-Fingerprint (pb_fingerprint/)
	1. Make directory in gb_id/t2avg_100/: mkdir pb_pkls/ 
	2. This code considers the atoms in the range of dis_width=7*latticxe_parameter above and below the dislocation position in y direction. dis_width" is defined in dir_names.py file.
	3. How to run: ovitos pb_fingerprint.py
III. Calculate descriptors for local atomic rearrangements (yprop_codes/)
	1. Make directory in gb_id/t2avg_100/: mkdir t2r_400/
		a. Make directory in gb_id/t2avg_100/t2r_400/: mkdir dmin
			a-1. Make directory in gb_id/t2avg_100/t2r_400/dmin/: mkdir dmin_pkls, mkdir dmin_dumps, mkdir dmin_svm
		b. Make directory in gb_id/t2avg_100/t2r_400/: mkdir phop
			b-1: Make directory in gb_id/t2avg_100/t2r_400/phop/: mkdir phop_pkls, mkdir phop_dumps, mkdir phop_svm
		c. Make directory in gb_id/t2avg_100/t2r_400/: mkdir gmean
			c-1: Make directory in gb_id/t2avg_100/t2r_400/gmean/: mkdir gmean_pkls, mkdir gmean_dumps, mkdir gmean_svm
	2. How to run: ovitos calc_dmin.py

Create training data-set
IV. CUR decomposition (cur_decomposition)
	1. Download pymf from https://github.com/cthurau/pymf
	2. Training data is written in yprop/yprop_svm/svm_train.pkl
	3. How to run: python2 yprop_cur.py
V. SVM Training




