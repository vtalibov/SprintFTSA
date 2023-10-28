# Vladimir O. Talibov
# Script for bulk processing of thermal shift assay data and thermal
# denaturation curves in general.
# Uses simple GUI to prompt input and output files.

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd
from os import path
import tkinter as tk
from tkinter import filedialog

def five_parametric_logistic_fixed(temp, infl, hill, assym):
    return min + ((max - min) / (1 + np.exp(hill*(infl - temp)))**assym)

# def four_parametric_logistic_fixed(temp, tm, hill):
#    return min + ((max - min)/ (1 + np.exp(hill*(tm - temp))))

def get_tm(infl, hill, assym):
    # Quantifies dTm value (midpoint transition) from the infliction
    # point for 5-parametric model.
    return infl - np.log(2**(1/assym) - 1)/ hill

def constrains(fluorescence):
    # Identifies transition region
    min = np.min(fluorescence[:np.argmax(np.diff(fluorescence))+1])
    min_ind = np.argmin(fluorescence[:np.argmax(np.diff(fluorescence))+1])
    max = np.max(fluorescence[np.argmax(np.diff(fluorescence)):])
    max_ind = np.argmax(fluorescence[np.argmax(np.diff(fluorescence)):]) + np.argmax(np.diff(fluorescence))
    return min, min_ind, max, max_ind

def normalize_fluorescence(fluorescence):
    min, min_ind, max, max_ind = constrains(fluorescence)
    return [((value - min)/(max - min))*100 for value in fluorescence]

def make_plot(temp, fluorescence, popt1, min, max):
    plt.xlabel("Temperature ($\degree$C)")
    plt.ylabel("Fluorescence")
    plt.title(column)
    plt.plot(temp, five_parametric_logistic_fixed(temp, *popt1), "--", color = 'black', label = "5PL")
    plt.vlines(get_tm(*popt1),min,max,color = "grey", 
               label = "$T_m$ = {} (5PL)".format(round(get_tm(*popt1),1)))
    plt.legend()
    plt.scatter(temp, fluorescence, s = 5)
    # to save plot to the same directory as the results file
    plot_path = path.join(path.dirname(output_filepath), "{}.png".format(column))
    plt.savefig(plot_path)
    plt.clf()

if __name__ == '__main__':

    version = '2310'

    # Some basic GUI
    win=tk.Tk()

    win.geometry("600x200")

    tk.Label(win, text="SprintTSA", font='arial 16 bold').pack(pady=15)
    tk.Label(win, text="Choose raw data file", font='arial 14 ').pack(pady=15)

    # Function to prompt for input and output
    def get_file_path():
        global input_filepath, output_filepath
        input_filepath = filedialog.askopenfilename(title = "Open raw data file",
                                                    filetypes=(("CSV Files","*.csv"),))
        output_filepath = filedialog.asksaveasfilename(title = "Save results table",
                                                       filetypes=(("CSV Files","*.csv"),))
        win.destroy() # not elegant, but to close tkinter window after retrieving path to the raw data.

    # Create a button to trigger the dialog
    button = tk.Button(win, text="Open", command=get_file_path)
    button.pack()

    tk.Label(win, text="{}".format(version), font='arial 12 ').pack(pady=15)

    win.mainloop()

    df = pd.read_csv(input_filepath)
    temp = df.iloc[:, 0]
    output=pd.DataFrame()
    for column in df.iloc[:,1:]:
        data=pd.DataFrame()
    #   Comment/uncomment lines below to fit&plot data with normalization.
    #    fluorescence = normalize_fluorescence(df[column]) # normalize data
        fluorescence = df[column] # do not normalize data
        min, min_ind, max, max_ind = constrains(fluorescence)
    #   Constrains for hill and assymetry coefficients are [-3, 3]; infliction point is
    #   fitted in a region between start and end of the transition.
        popt1, _ = opt.curve_fit(five_parametric_logistic_fixed, temp[min_ind : max_ind + 1], 
                               fluorescence[min_ind : max_ind +1], 
                               bounds = ([temp[min_ind],-3,-3],[temp[max_ind],3,3]))
        data['Well'] = [column]
        data['5PL'] = [get_tm(*popt1)]
        output = pd.concat([output, data])
        make_plot(temp, fluorescence, popt1, min, max)
    output.to_csv(output_filepath, index = False)
