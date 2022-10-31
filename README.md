# Calcium_ExTaSy
A public repository for implementing code for analyzing spontaneous calcium transients
suite2p installation w/in conda environment:

*NOTE: updates to suite2p have created some issues within the GUI making it less user friendly.
	This installation will allow for users to run suite2p and prevents errors from arrising
	IF someone is planning to do cell detection with suite2p, other errors might arrise from this
	installation method.

download suite2p-revert-435 branch .zip folder and move to appropriate location
Instead, you could also go directly to their Github and download suite2p v0.8.0, instead; either will produce usable results

in anaconda python terminal

>> cd "C:\users\*username*\suite2p-revert"

>> conda env create -f environment.yml
		'$' indicate code snippet starts

	This environment can now be activated using $ conda activate suite2p
	deactivate this environment with $ conda deactivate
	
All pip packages can be installed simultaneously by separating the packages to install with spaces between each. 
If you specify a version, you CANNOT have any spaces between the name of the program and the version 
(e.g. python=3.7 works; python = 3.7 will fail; == is necessary for packages)

>> pip install BaselineRemoval 
	to do interative baseline subtraction on fluorescent data

>> conda install jupyter notebook matplotlib pandas numba h5py natsort scipy seaborn

>> pip install scanimage-tiff-reader tifffile tdqm

>> pip install PyQt5==5.15.1 PyQt5.sip==12.8.1 pyqtgraph==0.11.0
 ***These versions will prevent an error pertailing to PyQt5.QtGui / 'QDialog' from arising
	other versions of packages do not matter

pip install rastermap cellpose
	***rastermap is necessary and allows for visualization of activity over time graphically
		cellpose is developed by suite2p devs and designed to detect somas...it is therefore optional
		
open s2p GUI: $ python -m suite2p
