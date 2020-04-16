import os
import sys

print("Please check the three scripts are consistent and writing and reading from the same directory. Continue y/n?")
proceed = input()
if proceed != 'y':
	sys.exit()

os.system('python simulationPatchyProteinMS.py')
os.system('python discretizePatchyProtein.py')
os.system('source activate base')
#os.system('source activate pyemma')
os.system('python generateMSMpatchyProtein.py')
