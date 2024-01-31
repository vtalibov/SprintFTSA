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
from tkinter import ttk


def constrains(fluorescence):
    '''
    Identifies transition region, assuming the temperature gradient is linear.
    '''
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
    '''
    Quantifies dTm value (midpoint transition) from the infliction point for 
    5-parametric sigmoidal model.
    :infl: - infliction point of the curve
    :hill: - coefficient of the curve steepness
    :assym: - curve asymmetry coefficient. If 1, infliction point = midpoint.
    '''
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
            tm_analysis_button["state"] = "normal"
            isothermal_analysis_button["state"] = "normal"

        
        def open_data(input):
            data = np.genfromtxt(input, delimiter=',', skip_header=1)
            with open(input, 'r') as file:
                header = file.readline().strip().split(',')
            return data, header
        
        def tm_analysis():
            experimental_data, header = open_data(input_filepath)
            temperature = experimental_data[:, 0] # retrieve temperature values
            output = np.array(['Well', 'infl', 'hill', 'assym', 'Tm'])
            for well in range(1, experimental_data.shape[1]):
                if normalize_checkbox_value.get() == 1:
                    fluorescence = normalize_fluorescence(experimental_data[:, well])
                else:
                    fluorescence = experimental_data[:, well]
                fit_parameters = fit_curve(temperature, fluorescence)
                if fit_parameters is None:
                    results = np.array([header[well], np.nan, np.nan, np.nan, np.nan])
                else:
                    results = np.array([header[well], fit_parameters[0], fit_parameters[1],
                                       fit_parameters[2], get_tm(*fit_parameters)])
                    title = header[well]
                    make_plot(temperature, fluorescence, fit_parameters, min_fluor, max_fluor, title)
                output = np.vstack((output, results))
            np.savetxt(output_filepath, output, fmt="%s", delimiter=',')

#       "end-1c" to avoid new line.
        def retrieve_isothermal_temperature_input():
            return isothermal_temperature_input.get("1.0","end-1c")

        def isothermal_analysis():
            try:
                temperature_of_interest = float(retrieve_isothermal_temperature_input())
            except Exception:
                tk.messagebox.showerror("Error", "Invalid temperature value.")
            experimental_data, header = open_data(input_filepath)
            temperature = experimental_data[:, 0] # retrieve temperature values
            output = np.array(['Well', '%Unfolded at {}C'.format(temperature_of_interest)])
            for well in range(1, experimental_data.shape[1]):
                fluorescence = normalize_fluorescence(experimental_data[:, well])
                fit_parameters = fit_curve(temperature, fluorescence)
                if fit_parameters is None:
                    results = np.array([header[well], np.nan])
                else:
                    min_fluor = 0
                    max_fluor = 100
                    results = np.array([header[well], 
                                        five_parametric_logistic_fixed(temperature_of_interest,
                                                                       *fit_parameters)])
                    output = np.vstack((output, results))
            np.savetxt(path.join(path.dirname(output_filepath), "isothermal_analysis.csv"), 
                       output, fmt="%s", delimiter=",")
            
#       interface objects
        picker_button = tk.Button(master, text="Open", command=io_path)
        tm_analysis_button = tk.Button(master, text="Analyse", state='disabled',
                                       command=tm_analysis)
        isothermal_analysis_button = tk.Button(master, text="Analyse", state='disabled',
                                               command=isothermal_analysis)
        isothermal_temperature_input = tk.Text(master, height = 1, width = 10)
        exit_button = tk.Button(master, text='Exit', command=quit)
        
        normalize_checkbox_value = tk.IntVar()
        normalize_checkbox = tk.Checkbutton(master, text='Normalize data',
                                            variable=normalize_checkbox_value, onvalue=1, offvalue=0)
        grid_output_checkbox_value = tk.IntVar()
        grid_output_checkbox = tk.Checkbutton(master, text='Microtiter plate grid',
                                              variable=grid_output_checkbox_value, onvalue=1,
                                              offvalue=0)
        separator_1 = ttk.Separator(root, orient='horizontal')
        separator_2 = ttk.Separator(root, orient='horizontal')
        
#       setting title and grid
        root.title("SprintFTSA")
        root.resizable(width=True, height=True)
        root.columnconfigure(0, weight=1)
        root.columnconfigure(1, weight=3)

#       interface
        tk.Label(master, text="SprintFTSA", font='arial 16 bold').grid(row = 0, column = 0,
                                                                       columnspan=2, pady=15)
        tk.Label(master, text=version, font='arial 12').grid(row=1, column=0, columnspan=2)
        picker_button.grid(row=2, column=0, columnspan=2, pady=5)
        separator_1.grid(row=3, column=0, columnspan=2, sticky='we')
        tk.Label(master, text='Tm analysis', font='Arial 14').grid(row=3, column=0, columnspan=2,
                                                                   pady=10)
        normalize_checkbox.grid(row=4, column=0, sticky='e', pady=5)
        grid_output_checkbox.grid(row=4, column=1, pady=5)
        tm_analysis_button.grid(row=5, column=0, columnspan=2, pady=5)
        separator_2.grid(row=6, column=0, columnspan=2, sticky='we')
        tk.Label(master, text='Isothermal analysis', font='Arial 14').grid(row=6, column=0,
                                                                           columnspan=2, pady=10)
        tk.Label(master, text='Enter temperature', font='Arial 10').grid(row=7, column = 0,
                                                                         sticky='e', pady=5)
        isothermal_temperature_input.grid(row=7, column=1, pady=5)
        isothermal_analysis_button.grid(row=8, column=0, columnspan=2, pady=5)
        exit_button.grid(row=9, column=0, columnspan=2, pady=15)

version = 'v.2312'

if __name__ == '__main__':
    root = tk.Tk()
    app = Application(root)
    root.title('SprintFTSA')
    root.minsize(300,300)
    root.mainloop()