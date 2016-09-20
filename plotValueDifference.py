import sys
import numpy as np
import scipy
import plot as Plot
import glob
import os
import logging
from serializer import Serializer
import matplotlib.pyplot as plt
import sets
import random
import argparse
#from emd import emd
from tabulate import tabulate
import fileinput

class ValueClass:
  def __init__(self):
    pass
  
  numObstacles = 0
  
  processCovariance = 0.0
  
  observationCovariance = 0.1
  
  valueABT = 0.0
  
  valueLQG = 0.0
  
  def serialize(self, directory):
    filename = ("out_" + 
               str(self.processCovariance) + 
               "_" + str(self.observationCovariance) + 
               ".txt")    
    files = glob.glob(os.path.join(directory + "/" + str(self.numObstacles) + "/" + filename))    
    if len(files) > 0: 
      data = 0
      with open(files[0], 'r') as f:
        # read a list of lines into data
        data = f.readlines()      
      for i in xrange(len(data)):	
        if "numObstacles" in data[i]:
	  data[i] = "numObstacles = " + str(self.numObstacles) + " \n"
	elif "processCovariance" in data[i]:
	  data[i] = "processCovariance = " + str(self.processCovariance)+ " \n"
	elif "observationCovariance" in data[i]:
	  data[i] = "observationCovariance = " + str(self.observationCovariance)+ " \n"
	elif "valueABT" in data[i] and self.valueABT != 0:
	  data[i] = "valueABT = " + str(self.valueABT) + " \n"
	elif "valueLQG" in data[i] and self.valueLQG != 0:
	  data[i] = "valueLQG = " + str(self.valueLQG) + " \n"
      with open(files[0], 'w') as f:
        f.writelines(data)      
    else:
      path = directory + "/" + str(self.numObstacles) + "/" + filename      
      with open(path, 'a+') as f:
	f.write("numObstacles = " + str(self.numObstacles) + " \n \n")
	f.write("processCovariance = " + str(self.processCovariance) + " \n \n")
	f.write("observationCovariance = " + str(self.observationCovariance) + " \n \n")
	f.write("valueABT = " + str(self.valueABT) + " \n \n")
	f.write("valueLQG = " + str(self.valueLQG) + " \n \n")

class PlotStats:
    def __init__(self, directory, loadValuesFromFiles):
        if directory == None:
	  print "No directory given"
	  return
        valueClasses = self.getValues(directory, loadValuesFromFiles)
        for v in valueClasses:
	  v.serialize(directory)
        self.generatePlots(valueClasses)
        sleep
        
    def generatePlots(self, valueClasses):
      motionErrors = []
      for valueClass in valueClasses:
	if not valueClass.processCovariance in motionErrors:
	  motionErrors.append(valueClass.processCovariance)
      for motionError in motionErrors:
	'''
	Generate a plot for each motion error
	'''
	valueClassesWithMotionError = [v for v in valueClasses 
				       if v.processCovariance == motionError]
	observationErrors = [v.observationCovariance for v in valueClassesWithMotionError]
	plotSets = []
	for observationError in observationErrors:
	  """
	  Generate a line for each observation error
	  """
	  valueClassesWithObservationError = [v for v in valueClassesWithMotionError
				              if v.observationCovariance == observationError]
	  setObservationError = []
	  for v in valueClassesWithObservationError:
	    absValue = 0
	    if v.valueABT != 0:
	      absValue = (v.valueABT - v.valueLQG) / np.fabs(v.valueABT)
	    setObservationError.append([v.numObstacles, absValue, 0.0])
	  plotSets.append(np.array(setObservationError))
	self.plotSets(plotSets, "dunnoyet.pdf", motionError, observationErrors)	
	
    def plotSets(self, sets, filename, motionError, observationErrors):
        labels = []
        for observationError in observationErrors:
	  labels.append(r'$err_{obs} = $' + str(observationError))
	plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=25)
	fig = plt.figure(figsize=(6, 4))
	ax = plt.subplot(1,1,1)
	for i in xrange(len(sets)):
	  ax.errorbar(sets[i][:,0], 
		      sets[i][:,1], 
		      yerr=sets[i][:,2],
		      label=labels[i], 
		      linewidth=5,
		      linestyle='solid')
	box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        
        plt.xlabel(r'\# obstacles')
        plt.ylabel(r'$\mathbf{\frac{V_{ABT} - V_{LQG-MP}}{\left |V_{ABT} \right |}}$', fontsize=50)
        plt.xlim([0, 30])
        plt.ylim([-1, 1])
        plt.xticks(np.linspace(0, 30, 4))
        plt.show()
        plt.clf()
		  
	sleep
	
    def getValues(self, directory, loadValuesFromFiles):       
        obstacleList = [0, 10, 20, 30]
        valueClasses = []
        if loadValuesFromFiles:
          valueClasses = self.loadValueClassesFromFiles(directory, obstacleList)
        for numObstacles in obstacleList:	  	  
	  obstacleDirectory = directory + "/" + str(numObstacles)
	  logFiles = glob.glob(os.path.join(obstacleDirectory + "/log_*"))
	  for logFile in logFiles:
	    secondLogFileParts = logFile.split("/")
	    secondLogFileParts = secondLogFileParts[-1].split("_")
	    processCovariance = float(secondLogFileParts[3])
	    observationCovarianceString = secondLogFileParts[4]
	    if "log" in secondLogFileParts[4]:
	       observationCovarianceString = secondLogFileParts[4].split(".log")[0]
	    observationCovariance = float(observationCovarianceString)
	    v = 0
	    valueClassExists, index = self.valueClassExists(valueClasses,
			                                    numObstacles,
			                                    processCovariance,
			                                    observationCovariance)
	    if valueClassExists:
	      v = valueClasses[index]
	    else:	      
	      v = ValueClass()
	      v.numObstacles = numObstacles
	      v.processCovariance = processCovariance
	      v.observationCovariance = observationCovariance
	     
	    meanReward = self.calcMeanReward(logFile)	    
	    logFileParts = logFile.split("/")
	    if "lqg" in logFileParts[-1]:
              v.valueLQG = meanReward
	    elif "abt" in part:
	      v.valueABT = meanReward
	    if not valueClassExists:
	      valueClasses.append(v)
	return valueClasses
      
    def loadValueClassesFromFiles(self, directory, obstacleList):      
      valueClasses = []
      for numObstacles in obstacleList:
	files = glob.glob(os.path.join(directory + "/" + str(numObstacles) + "/out_*"))
	data = 0
	for file in files:
	  v = ValueClass()
	  with open(file, 'r') as f:
	    data = f.readlines()
	  for i in xrange(len(data)):
	    if "numObstacles" in data[i]:
	      v.numObstacles = float(data[i].strip().split(" = ")[1])
	    elif "processCovariance" in data[i]:
	      v.processCovariance = float(data[i].strip().split(" = ")[1])
	    elif "observationCovariance" in data[i]:
	      v.observationCovariance = float(data[i].strip().split(" = ")[1])
	    elif "valueABT" in data[i]:
	      v.valueABT = float(data[i].strip().split(" = ")[1])
	    elif "valueLQG" in data[i]:
	      v.valueLQG = float(data[i].strip().split(" = ")[1])
	  print v.numObstacles
	  print v.processCovariance
	  print v.observationCovariance
	  print v.valueABT
	  print v.valueLQG; sleep
	  valueClasses.append(v)
	    
	
    def valueClassExists(self, 
			 valueClasses,
			 numObstacles,
			 processCovariance, 
			 observationCovariance):      
      for i in xrange(len(valueClasses)):
	if (valueClasses[i].numObstacles == numObstacles and
	    valueClasses[i].processCovariance == processCovariance and 
            valueClasses[i].observationCovariance == observationCovariance):
	  return True, i
      return False, 0
	  
	 
    def calcMeanReward(self, file):
        print "calc mean reward..."
        reward = 0
        rewardsPerRun = []
        with open(file, "r") as f:
	  for line in f:
	    if "R: " in line:
	      reward += float(line.split("R: ")[-1].strip())	      
	    elif ("RUN #" in line or 
	          "Run #" in line or
		  "#############" in line):
	      rewardsPerRun.append(reward)
	      reward = 0
	sumRewards = sum(rewardsPerRun)	
	if len(rewardsPerRun) == 0:
	  return 0.0
	return sumRewards / len(rewardsPerRun)
	  
	
    
    def calcRelativeValueDifference(self, dir):
	outFiles = glob.glob(os.path.join(dir, "out_*.txt")) 
	print outFiles
	
    
    
    
if __name__ == "__main__":
    logging_level = logging.DEBUG
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging_level)
    parser = argparse.ArgumentParser(description='LQG-MP plotting')    
    parser.add_argument("-d", "--directory", help="Path to the robot model file") 
    parser.add_argument("-l", "--loadValuesFromFiles", 
			help="Load the values from generated files",
			action="store_true") 
    args = parser.parse_args()
    PlotStats(args.directory, args.loadValuesFromFiles)
