# LaxWendroff

Run.py - python file that actually runs the Serre equation solver, current example is the soliton.

Makefile - make file to produce the c python linked with swig to be able to interface with python. Will require a proper link to Python.h, for me this was accomplished with '-I/home/jp/anaconda2/include/python2.7', (https://stackoverflow.com/questions/25882150/python-h-not-found-using-swig-and-anaconda-python). May have to do different things depending on your installation.
