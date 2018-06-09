import numpy as np
import sympy as sym
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats
import os, sys
from math import *
from bokeh.layouts import gridplot
from bokeh.plotting import output_file, figure, show
from bokeh.palettes import inferno
from bokeh.models import Label
from bokeh.models import Range1d

path = '/home/yagafarova/DNA_packaging/dsDNA/'

#Langevin vs. Nose-Hoover
p = figure(x_axis_label='Radius of compression, nm', y_axis_label='Pressure, atm', 
           plot_width=800, plot_height=600,
           x_range = (9.6, 20.4), y_range = (-20, 520))
output_file('len_nh.html', title="len_nh")
p.xaxis.axis_label_text_font_size = "22pt"
p.xaxis.major_label_text_font_size = "18pt"
p.yaxis.axis_label_text_font_size = "22pt"
p.yaxis.major_label_text_font_size = "18pt"
p.legend.label_text_font_size = "18pt"

files = ['VelocityVerlet_Langevin/dsDNA_4000bp_mM50/push/0.1/dsDNA_4000bp_mM50_pressure_push.dat', 
         '4000bp_concentration/dsDNA_4000bp_mM50/push_10k_steps/dsDNA_4000bp_mM50_pressure_push.dat',
         '4000bp_concentration/dsDNA_4000bp_mM50/push_100k_steps/dsDNA_4000bp_mM50_pressure_push.dat']
labels = ['Langevin',
          'Nose-Hoover: 10 M steps',
          'Nose-Hoover: 100 M steps']
colors = inferno(len(labels)+1)
for fname,l,c in zip(files, labels, colors):
    fname = path + fname
    radius = [float(i.split('\t')[1].strip()) for i in open(fname)]
    pressure = [float(i.split('\t')[2].strip())*1.6*9.87 for i in open(fname)]
    #start = np.where(np.array(radius)<20)[0]
    #p.line(np.array(radius)[start], np.array(pressure)[start], color= c, legend=l, line_width=2)
    p.line(radius, pressure, color= c, legend=l, line_width=2)
show(p)