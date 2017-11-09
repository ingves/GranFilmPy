import os
import sys
import matplotlib.pyplot as plt

# --- Set essential paths

gf_package = os.path.dirname(os.path.realpath(__file__)) # Path to the GranFilm_Package directory
python_interface = os.path.join(gf_package,'Python_Interface') # Path to the directory of the .py files constituting the python interface for GranFilm 
input_file = os.path.join(gf_package,'input_file') # Path to the input file for GranFilm simulation

# --- Import the Python interface

sys.path.append(python_interface)
from GranFilmPy import GranFilm

gf = GranFilm()
#gf()
print("Classic : OK")

gf = GranFilm(points_file="classic")
#gf()
#gf["Potential"]()
print("Classic + potential : OK")

gf = GranFilm(points_file='classic',radius_ratios=[1.0,0.7,0.4],media="air,mgo,ag,ag,mgo,mgo,ag,ag")
#gf()
#gf["Potential"]()
print("coated + potential : OK")

gf = GranFilm(points_file='classic',radius=[8.0,5.0],theta=0.0)
gf()
gf["Potential"]()
print("oblate + potential : OK")

gf = GranFilm(points_file='classic',radius=[5.0,8.0],theta=0.0)
gf()
gf["Potential"]()
print("prolate+ potential : OK")
