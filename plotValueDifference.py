import sys
import numpy as np
import scipy
import plot as Plot
import glob
import os
import logging
from serializer import Serializer
#from scipy.stats import multivariate_normal
#from librobot import *
#from libutil import *
#from EMD import EMD
import sets
import random
import argparse
#from emd import emd
from tabulate import tabulate

class PlotStats:
    def __init__(self, dir, save_plots, show_particles, plot_emds, collision_is_failed):
	self.calcRelativeValueDifference(dir)
    
    def calcRelativeValueDifference(self, dir):
	outFiles = glob.glob(os.path.join(dir, "out_*.txt")) 
	print outFiles
	
    
    
    
if __name__ == "__main__":
    logging_level = logging.DEBUG
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging_level)
    parser = argparse.ArgumentParser(description='LQG-MP plotting')    
    parser.add_argument("-d", "--directory", help="Path to the robot model file")
    parser.add_argument("-s", "--save", 
                        help="Save the plots", 
                        action="store_true")
    parser.add_argument("-p", "--particles", 
                        help="Show particles", 
                        action="store_true")
    parser.add_argument("-e", "--emd", 
                        help="Plot emds", 
                        action="store_true")
    parser.add_argument("-cf", "--collision_is_failed",
                        help="A run in which the robot collides, counts as an unsuccessful run",
                        action="store_true")
    args = parser.parse_args() 
    PlotStats(args.directory, args.save, args.particles, args.emd, args.collision_is_failed)