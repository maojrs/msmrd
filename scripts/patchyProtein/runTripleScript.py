import os
import sys

print("Please check the three scripts are consistent and writing and reading from the same directory. To continue press y")
proceed = input()
if proceed != 'y':
	sys.exit()

os.system('python simulationPatchyProteinMS.py')
os.system('python discretizePatchyProtein.py')
os.system('python generateMSMpatchyProtein.py')
