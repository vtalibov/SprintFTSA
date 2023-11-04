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


def constrains(fluorescence):
    # Identifies transition region
    min_fluor = np.min(fluorescence[:np.argmax(np.diff(fluorescence))+1])
    min_ind = np.argmin(fluorescence[:np.argmax(np.diff(fluorescence))+1])
    max_fluor = np.max(fluorescence[np.argmax(np.diff(fluorescence)):])
    max_ind = np.argmax(fluorescence[np.argmax(np.diff(fluorescence)):]) + np.argmax(np.diff(fluorescence))
    return min_fluor, min_ind, max_fluor, max_ind

def normalize_fluorescence(fluorescence):
    min_fluor, min_ind, max_fluor, max_ind = constrains(fluorescence)
    return [((value - min_fluor)/(max_fluor - min_fluor))*100 for value in fluorescence]

def five_parametric_logistic_fixed(temp, infl, hill, assym):
    return min_fluor + ((max_fluor - min_fluor) / (1 + np.exp(hill*(infl - temp)))**assym)

def get_tm(infl, hill, assym):
    # Quantifies dTm value (midpoint transition) from the infliction
    # point for 5-parametric model.
    return infl - np.log(2**(1/assym) - 1)/ hill

def make_plot(temp, fluorescence, popt1, min, max):
    plt.xlabel("Temperature ($\degree$C)")
    plt.ylabel("Fluorescence")
    plt.title(well_name)
    plt.plot(temp, five_parametric_logistic_fixed(temp, *popt1), "--", color = 'black', label = "5PL")
    plt.vlines(get_tm(*popt1),min_fluor,max_fluor,color = "grey", 
               label = "$T_m$ = {} (5PL)".format(round(get_tm(*popt1),1)))
    plt.legend()
    plt.scatter(temp, fluorescence, s = 5)
    # to save plot to the same directory as the results file
    plot_path = path.join(path.dirname(output_filepath), "{}.png".format(well_name))
    plt.savefig(plot_path)
    plt.clf()    

# Main class
class Application():
    def __init__(self, master):
        def io_path():
            global input_filepath, output_filepath
            input_filepath = filedialog.askopenfilename(title = "Open raw data file",
                                                        filetypes=(("CSV Files","*.csv"),))
            output_filepath = filedialog.asksaveasfilename(title = "Save results table",
                                                        filetypes=(("CSV Files","*.csv"),))
            analyse_button["state"] = "normal"
        def analyse():
            # non-elegant to pass non-declared variables to all functions
            global well_name, min_fluor, min_ind, max_fluor, max_in
            df = pd.read_csv(input_filepath)
            temp = df['Temperature']
            output=pd.DataFrame()
            for well_name in df.iloc[:,1:]:
                data=pd.DataFrame()
                fluorescence = df[well_name]
                min_fluor, min_ind, max_fluor, max_ind = constrains(fluorescence)
            # Constrains for hill and assymetry coefficients are [-3, 3]; infliction point is
            # fitted in a region between start and end of the transition.
                popt1, _ = opt.curve_fit(five_parametric_logistic_fixed, temp[min_ind : max_ind + 1], 
                                         fluorescence[min_ind : max_ind +1], 
                                         bounds = ([temp[min_ind],-3,-3],[temp[max_ind],3,3]))
                data['Well'] = [well_name]
                data['5PL'] = [get_tm(*popt1)]
                output = pd.concat([output, data])
                make_plot(temp, fluorescence, popt1, min_fluor, max_fluor)
            output.to_csv(output_filepath, index = False)      
        
        frame = tk.Frame(master)
        frame.pack()

        picker_button = tk.Button(master, text="Open", command=io_path)
        analyse_button = tk.Button(master, text="Analyse", state='disabled', command=analyse)
        exit_button = tk.Button(master, text='Exit', command=quit)
        
        normalize_checkbox_value = tk.BooleanVar()
        normalize_checkbox = tk.Checkbutton(master, text='Normalize data',
                                            variable=normalize_checkbox_value,
                                            onvalue=True, offvalue=False)
        
        tk.Label(master, text="SprintFTSA", font='arial 16 bold').pack(pady=15)
        tk.Label(master, text=version, font='arial 12').pack()

        normalize_checkbox.pack(pady=15)
        picker_button.pack(pady=15)
        analyse_button.pack(pady=15)
        exit_button.pack(side=tk.BOTTOM)

version = '2311'

if __name__ == '__main__':
    root = tk.Tk()
    app = Application(root)
    root.title('SprintFTSA')
    root.minsize(300,300)
    root.mainloop()