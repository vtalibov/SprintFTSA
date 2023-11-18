# Vladimir O. Talibov
# Script for bulk processing of thermal shift assay data and thermal
# denaturation curves in general.
# Uses simple GUI to prompt input and output files.

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
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

def five_parametric_logistic_fixed(temperature, infl, hill, assym):
    return min_fluor + ((max_fluor - min_fluor) / (1 + np.exp(hill*(infl - temperature)))**assym)

def fit_curve(temperature, fluorescence):
    global min_fluor, min_ind, max_fluor, max_ind
    min_fluor, min_ind, max_fluor, max_ind = constrains(fluorescence)
    try:
        # Constrains for hill and assymetry coefficients are [-3, 3]; infliction point is
        # fitted in a region between start and end of the transition.
        parameters, covariance = opt.curve_fit(five_parametric_logistic_fixed, 
                                               temperature[min_ind : max_ind + 1],
                                               fluorescence[min_ind : max_ind +1],
                                               bounds = ([temperature[min_ind],-3,-3],
                                                         [temperature[max_ind],3,3]))
        if np.any(np.diag(covariance) <= 0): # quality of the fit check
            raise RuntimeError("No convergence.")
        return parameters
    except Exception as error:
        print("Curve fitting error {}".format(error))
        return None

def get_tm(infl, hill, assym):
    # Quantifies dTm value (midpoint transition) from the infliction
    # point for 5-parametric model.
    return infl - np.log(2**(1/assym) - 1)/ hill

def make_plot(temperature, fluorescence, fit_parameters, min_fluor, max_fluor, plot_title = 'Title'):
    plt.xlabel("Temperature ($\degree$C)")
    plt.ylabel("Fluorescence")
    plt.title(plot_title)
    plt.plot(temperature, five_parametric_logistic_fixed(temperature, *fit_parameters),
             "--", color = 'black', label = "5PL")
    plt.vlines(get_tm(*fit_parameters), min_fluor, max_fluor, color = "grey",
               label = "$T_m$ = {} (5PL)".format(round(get_tm(*fit_parameters),1)))
    plt.legend()
    plt.scatter(temperature, fluorescence, s = 5)
    # to save plot to the same directory as the results file
    plot_path = path.join(path.dirname(output_filepath), "{}.png".format(plot_title))
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
                                                           defaultextension=".csv",
                                                           filetypes=(("CSV Files","*.csv"),))
            analyse_button["state"] = "normal"
        def analyse():
            experimental_data = np.genfromtxt(input_filepath, delimiter=',', skip_header=1)
            with open(input_filepath, 'r') as file:
                header = file.readline().strip().split(',')
            temperature = experimental_data[:, 0] # retrieve temperature values
            output = np.array(['Well', '5PL'])
            for well in range(1, experimental_data.shape[1]):
                if normalize_checkbox_value.get() == 1:
                    fluorescence = normalize_fluorescence(experimental_data[:, well])
                else:
                    fluorescence = experimental_data[:, well]
                fit_parameters = fit_curve(temperature, fluorescence)
                if fit_parameters is None:
                    results = np.array([header[well], np.nan])
                else:
                    results = np.array([header[well], get_tm(*fit_parameters)])
                    title = header[well]
                    make_plot(temperature, fluorescence, fit_parameters, min_fluor, max_fluor, title)
                output = np.vstack((output, results))
            np.savetxt(output_filepath, output, fmt="%s", delimiter=',')

        frame = tk.Frame(master)
        frame.pack()

        picker_button = tk.Button(master, text="Open", command=io_path)
        analyse_button = tk.Button(master, text="Analyse", state='disabled', command=analyse)
        exit_button = tk.Button(master, text='Exit', command=quit)
        
        normalize_checkbox_value = tk.IntVar()
        normalize_checkbox = tk.Checkbutton(master, text='Normalize data',
                                            variable=normalize_checkbox_value, onvalue=1, offvalue=0)
        
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