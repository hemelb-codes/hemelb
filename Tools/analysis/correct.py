import os



os.system("python stress.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715131802/")

os.system("python allPlanes.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132215/")
os.system("python stress.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132215/")
os.system("python velocityMag.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132215/")
os.system("python velocityField.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132215/")


os.system("python allPlanes.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132629/")
os.system("python stress.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132629/")
os.system("python velocityMag.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132629/")
os.system("python velocityField.py ~/devel/ccs/fabric/hemelb/results/nhnn_jun2014_archer_184501966cbb+_20140715132629/")

