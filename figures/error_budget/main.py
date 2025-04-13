import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import cmath

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

def make_autopct(allvals):
  def autopct(pct):
    total = sum(allvals)
    val = pct * total / 100.0
    # Customize the returned string as desired
    return f"${val:.2f}\%$"
  return autopct

def main():
  with open('../style.py', 'r') as f:
    s = f.read()
    plot = eval(s)

  # Data
  labels = ["$\mathrm{Scale}$", "$\mathrm{Electroweak}$", "$\mathrm{PDF-MHOU}$", "$\mathrm{Top}$", "$\mathrm{Top-Bottom}$"]
  values = [1.95, 0.3, 0.3, 0.2, 1.42]

  colors = plot["color"][0:len(values)]

  # Create the pie chart
  plt.figure(figsize=(4,4))
  plt.pie(values, labels=labels, autopct=make_autopct(values), startangle=140, colors=colors, pctdistance=0.8, textprops={'color':"w"}, wedgeprops = { 'linewidth' : 0.8, 'edgecolor' : 'w' , 'joinstyle' : 'round'})

  # Add a title and make the plot look a bit nicer
  plt.axis('equal')  # Ensures the pie chart is circular
  plt.tight_layout()
  plt.savefig("../../Images/error_budget_before.pdf")

  # Show the pie chart
  plt.show()

  # Data
  labels = ["$\mathrm{Scale}$", "$\mathrm{Electroweak}$", "$\mathrm{PDF-MHOU}$", "$\mathrm{Top}$", "$\mathrm{Top-Bottom}$"]
  values = [1.95, 0.3, 0.3, 0.2, 0.2]

  colors = plot["color"][0:len(values)]

  # Create the pie chart
  plt.figure(figsize=(4,4))
  plt.pie(values, labels=labels, autopct=make_autopct(values), startangle=140, colors=colors, pctdistance=0.8, textprops={'color':"w"}, wedgeprops = { 'linewidth' : 0.8, 'edgecolor' : 'w' , 'joinstyle' : 'round'})

  # Add a title and make the plot look a bit nicer
  plt.axis('equal')  # Ensures the pie chart is circular
  plt.tight_layout()
  plt.savefig("../../Images/error_budget_after.pdf")

  # Show the pie chart
  plt.show()

if __name__=="__main__":
  main()