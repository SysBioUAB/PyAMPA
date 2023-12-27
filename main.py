import tkinter as tk
import subprocess
import math
import random
import decimal
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import os
import sys
import tempfile
from tkinter import filedialog, messagebox
from tkinter import ttk
import tkinter.font as font
from tkinter.ttk import Progressbar
from PIL import ImageTk, Image
from Bio import SeqIO
#from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as BioIP
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.ensemble import RandomForestClassifier
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.backends.backend_agg as agg
import matplotlib.pylab as pylab
import matplotlib.patches as patches
import matplotlib.lines as lines
from contextlib import contextmanager

POPULATION_SIZE = 100
NUM_GENERATIONS = 100
MAX_NO_IMPROVEMENT = 20
CROSSOVER_RATE = 0.8
MUTATION_RATE = 0.2
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
amino_acids_mutate = ['A', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Define your weights
w1 = 1
w2 = 1
w3 = 1

# Load AMPValidate vectorizer
with open('AMPValidate.pkl', 'rb') as f:
    model_amp = pickle.load(f)
with open('amp_validate_vectorizer.pkl', 'rb') as f:
    vectorizer_amp = pickle.load(f)

# Load the Hemolysis model and vectorizer
with open("hemolysis_model.pkl", 'rb') as file:
    model_hemo = pickle.load(file)
with open("hemolysis_vectorizer.pkl", 'rb') as file:
    vectorizer_hemo = pickle.load(file)

def center_window(window):
    window.update_idletasks()
    width = window.winfo_width()
    height = window.winfo_height()
    x = (window.winfo_screenwidth() // 2) - (width // 2)
    y = (window.winfo_screenheight() // 2) - (height // 2)
    window.geometry('{}x{}+{}+{}'.format(width, height, x, y))

# Function to split peptide sequences into all possible subsequences of length 2
def split_sequence(sequence):
    return ' '.join([sequence[i:i+2] for i in range(len(sequence) - 1)])

def round_to_sig_figs(x, sig_figs):
    if x != 0:
        magnitude = int(math.floor(math.log10(abs(x))))
        if magnitude >= 0:
            # Number is greater or equal to 1
            decimal_places = sig_figs - magnitude - 1
            if decimal_places <= 0:
                return "{:.0f}".format(round(x))
            else:
                return "{:0.{}f}".format(round(x, decimal_places), decimal_places)
        else:
            # Number is less than 1
            return "{:0.{}f}".format(round(x, sig_figs - magnitude - 1), sig_figs - magnitude - 1)
    else:
        return "0"  # return "0" if x is 0


def calc_features(sequence):
    # Define the set of nonpolar residues
    nonpolar_residues = set('ACGILMFPWYV')

    # Calculate the percentage of nonpolar residues
    NP = sum(1 for aa in sequence if aa in nonpolar_residues) / len(sequence)

    # Calculate the presence of Trp and Tyr
    W = 1 if 'W' in sequence else 0
    Y = 1 if sequence.count('Y') >= 2 else 0

    # Calculate the isoelectric point
    protein_analysis = ProteinAnalysis(sequence)
    IP_val = protein_analysis.isoelectric_point()
    IP = 1 if IP_val > 10 else 0

    return NP, W, Y, IP

def calc_half_life(sequence):
    # Calculate the peptide features
    NP, W, Y, IP = calc_features(sequence)

    # Calculate the natural log of the half-life
    ln_t_half = 2.226 + (0.053 * NP * 100) - (1.515 * W) + (1.290 * Y) - (1.052 * IP)

    # Convert to the half-life
    t_half = math.exp(ln_t_half)

    return t_half


def fitness(antimicrobial_proba, hemolytic_proba, half_life, w1, w2, w3):
    # Convert half_life to a NumPy array
    half_life = np.array(half_life)

    # Limit the half-life to a maximum of 360 minutes
    half_life_limited = np.where(half_life > 360, 360, half_life)

    # Normalize the half-life to the range 0-1
    half_life_normalized = half_life_limited / 360

    # Calculate the fitness score as a weighted sum of the three factors
    return w1 * antimicrobial_proba - w2 * hemolytic_proba + w3 * half_life_normalized

def mutate_sequence(sequence):
    # Choose a random index for the mutation
    mutation_index = random.randint(0, len(sequence) - 1)

    # Choose a new amino acid that is different from the current one
    new_amino_acid = random.choice([aa for aa in amino_acids_mutate if aa != sequence[mutation_index]])

    # Create the new sequence with the mutation
    new_sequence = sequence[:mutation_index] + new_amino_acid + sequence[mutation_index + 1:]

    return new_sequence

# Probability of a sequence to be antimicrobial
def predict_proba_amp(sequence):

    # Split the sequence into all possible subsequences of length 2
    sequence_split = split_sequence(sequence)
    
    # Count the occurrences of each subsequence
    X = vectorizer_amp.transform([sequence_split])
    
    # Predict the probability using the trained model
    prediction_proba = model_amp.predict_proba(X)
    
    # Return the probability of the class being 1
    return prediction_proba[0][1]

# Probability of a sequence to be hemolytic
def predict_proba_hemo(sequence):
    
    # Split the sequence into all possible subsequences of length 2
    sequence_split = split_sequence(sequence)
    
    # Count the occurrences of each subsequence
    X = vectorizer_hemo.transform([sequence_split])
    
    # Predict the probability using the trained model
    prediction_proba = model_hemo.predict_proba(X)
    
    # Return the probability of the class being 1
    return prediction_proba[0][1]

def load_scale(scalename):
    """Method to load scale values for a given amino acid scale

    :param scalename: amino acid scale name, for available scales see the
        :class:`modlamp.descriptors.PeptideDescriptor()` documentation.
    :return: amino acid scale values in dictionary format.
    """
    # predefined amino acid scales dictionary
    scales = {
        'eisenberg': {'I': [1.4], 'F': [1.2], 'V': [1.1], 'L': [1.1], 'W': [0.81], 'M': [0.64], 'A': [0.62],
                      'G': [0.48], 'C': [0.29], 'Y': [0.26], 'P': [0.12], 'T': [-0.05], 'S': [-0.18], 'H': [-0.4],
                      'E': [-0.74], 'N': [-0.78], 'Q': [-0.85], 'D': [-0.9], 'K': [-1.5], 'R': [-2.5]}
    }

    if scalename in scales:
        return scales[scalename]
    else:
        raise ValueError(f"The scale {scalename} is not defined.")

def calculate_hydrophobic_moment(sequence, window=1000, angle=100):
    # if sequence is shorter than window, take the whole sequence instead
    window = min(window, len(sequence))
    
    # calculate descriptor values for each amino acid in the sequence
    d_eisberg = load_scale('eisenberg')
    descriptors = [d_eisberg[aa][0] for aa in sequence]
    
    # calculate the hydrophobic moment for each window in the sequence
    moments = []
    for i in range(len(descriptors) - window + 1):
        window_descriptors = descriptors[i:i + window]
        
        # calculate actual moment (radial)
        rads = angle * (np.pi / 180) * np.asarray(range(window))
        vcos = sum(window_descriptors * np.cos(rads))
        vsin = sum(window_descriptors * np.sin(rads))
        
        moment = np.sqrt(vsin**2 + vcos**2) / window
        moments.append(moment)
    
    # return the maximum hydrophobic moment
    return max(moments)    


def helical_wheel(sequence, colorcoding='rainbow', lineweights=True, filename=None, seq=False, moment=False):
    """A function to project a given peptide sequence onto a helical wheel plot. It can be useful to illustrate the
    properties of alpha-helices, like positioning of charged and hydrophobic residues along the sequence.

    :param sequence: {str} the peptide sequence for which the helical wheel should be drawn
    :param colorcoding: {str} the color coding to be used, available: *rainbow*, *charge*, *polar*, *simple*,
        *amphipathic*, *none*
    :param lineweights: {boolean} defines whether connection lines decrease in thickness along the sequence
    :param filename: {str} filename  where to safe the plot. *default = None* --> show the plot
    :param seq: {bool} whether the amino acid sequence should be plotted as a title
    :param moment: {bool} whether the Eisenberg hydrophobic moment should be calculated and plotted
    :return: a helical wheel projection plot of the given sequence (interactively or in **filename**)
    :Example:

    >>> helical_wheel('GLFDIVKKVVGALG')
    >>> helical_wheel('KLLKLLKKLLKLLK', colorcoding='charge')
    >>> helical_wheel('AKLWLKAGRGFGRG', colorcoding='none', lineweights=False)
    >>> helical_wheel('ACDEFGHIKLMNPQRSTVWY')

    .. image:: ../docs/static/wheel1.png
        :height: 300px
    .. image:: ../docs/static/wheel2.png
        :height: 300px
    .. image:: ../docs/static/wheel3.png
        :height: 300px
    .. image:: ../docs/static/wheel4.png
        :height: 300px

    .. versionadded:: v2.1.5
    """
    # color mappings
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    f_rainbow = ['#3e3e28', '#ffcc33', '#b30047', '#b30047', '#ffcc33', '#3e3e28', '#80d4ff', '#ffcc33', '#0047b3',
                 '#ffcc33', '#ffcc33', '#b366ff', '#29a329', '#b366ff', '#0047b3', '#ff66cc', '#ff66cc', '#ffcc33',
                 '#ffcc33', '#ffcc33']
    f_charge = ['#000000', '#000000', '#ff4d94', '#ff4d94', '#000000', '#000000', '#80d4ff', '#000000', '#80d4ff',
                '#000000', '#000000', '#000000', '#000000', '#000000', '#80d4ff', '#000000', '#000000', '#000000',
                '#000000', '#000000']
    f_polar = ['#000000', '#000000', '#80d4ff', '#80d4ff', '#000000', '#000000', '#80d4ff', '#000000', '#80d4ff',
               '#000000', '#000000', '#80d4ff', '#000000', '#80d4ff', '#80d4ff', '#80d4ff', '#80d4ff', '#000000',
               '#000000', '#000000']
    f_simple = ['#ffcc33', '#ffcc33', '#0047b3', '#0047b3', '#ffcc33', '#7f7f7f', '#0047b3', '#ffcc33', '#0047b3',
                '#ffcc33', '#ffcc33', '#0047b3', '#ffcc33', '#0047b3', '#0047b3', '#0047b3', '#0047b3', '#ffcc33',
                '#ffcc33', '#ffcc33']
    f_none = ['#ffffff'] * 20
    f_amphi = ['#ffcc33', '#29a329', '#b30047', '#b30047', '#f79318', '#80d4ff', '#0047b3', '#ffcc33', '#0047b3',
               '#ffcc33', '#ffcc33', '#80d4ff', '#29a329', '#80d4ff', '#0047b3', '#80d4ff', '#80d4ff', '#ffcc33',
               '#f79318', '#f79318']
    t_rainbow = ['w', 'k', 'w', 'w', 'k', 'w', 'k', 'k', 'w', 'k', 'k', 'k', 'k', 'k', 'w', 'k', 'k', 'k', 'k', 'k']
    t_charge = ['w', 'w', 'k', 'k', 'w', 'w', 'k', 'w', 'k', 'w', 'w', 'w', 'w', 'w', 'k', 'w', 'w', 'w', 'w', 'w']
    t_polar = ['w', 'w', 'k', 'k', 'w', 'w', 'k', 'w', 'k', 'w', 'w', 'k', 'w', 'k', 'k', 'k', 'k', 'w', 'w', 'w']
    t_simple = ['k', 'k', 'w', 'w', 'k', 'w', 'w', 'k', 'w', 'k', 'k', 'k', 'k', 'w', 'w', 'w', 'w', 'k', 'k', 'k']
    t_none = ['k'] * 20
    t_amphi = ['k', 'k', 'w', 'w', 'w', 'k', 'w', 'k', 'w', 'k', 'k', 'k', 'w', 'k', 'w', 'k', 'k', 'k', 'w', 'w']
    d_eisberg = load_scale('eisenberg')  # eisenberg hydrophobicity values for HM
    
    if lineweights:
        lw = np.arange(0.1, 5.5, 5. / (len(sequence) - 1))  # line thickness array
        lw = lw[::-1]  # inverse order
    else:
        lw = [2.] * (len(sequence) - 1)
    
    # check which color coding to use
    if colorcoding == 'rainbow':
        df = dict(zip(aa, f_rainbow))
        dt = dict(zip(aa, t_rainbow))
    elif colorcoding == 'charge':
        df = dict(zip(aa, f_charge))
        dt = dict(zip(aa, t_charge))
    elif colorcoding == 'polar':
        df = dict(zip(aa, f_polar))
        dt = dict(zip(aa, t_polar))
    elif colorcoding == 'simple':
        df = dict(zip(aa, f_simple))
        dt = dict(zip(aa, t_simple))
    elif colorcoding == 'none':
        df = dict(zip(aa, f_none))
        dt = dict(zip(aa, t_none))
    elif colorcoding == 'amphipathic':
        df = dict(zip(aa, f_amphi))
        dt = dict(zip(aa, t_amphi))
    else:
        print("Unknown color coding, 'rainbow' used instead")
        df = dict(zip(aa, f_rainbow))
        dt = dict(zip(aa, t_rainbow))
    
    # degree to radian
    deg = np.arange(float(len(sequence))) * -100.
    deg = [d + 90. for d in deg]  # start at 270 degree in unit circle (on top)
    rad = np.radians(deg)
    
    # dict for coordinates and eisenberg values
    d_hydro = dict(zip(rad, [0.] * len(rad)))
    
    # create figure
    fig = plt.figure(frameon=False, figsize=(3, 3))
    ax = fig.add_subplot(111)
    old = None
    hm = list()
    
    # iterate over sequence
    for i, r in enumerate(rad):
        new = (np.cos(r), np.sin(r))  # new AA coordinates
        if i < 18:
            # plot the connecting lines
            if old is not None:
                line = lines.Line2D((old[0], new[0]), (old[1], new[1]), transform=ax.transData, color='k',
                                    linewidth=lw[i - 1])
                line.set_zorder(1)  # 1 = level behind circles
                ax.add_line(line)
        elif 17 < i < 36:
            line = lines.Line2D((old[0], new[0]), (old[1], new[1]), transform=ax.transData, color='k',
                                linewidth=lw[i - 1])
            line.set_zorder(1)  # 1 = level behind circles
            ax.add_line(line)
            new = (np.cos(r) * 1.2, np.sin(r) * 1.2)
        elif i == 36:
            line = lines.Line2D((old[0], new[0]), (old[1], new[1]), transform=ax.transData, color='k',
                                linewidth=lw[i - 1])
            line.set_zorder(1)  # 1 = level behind circles
            ax.add_line(line)
            new = (np.cos(r) * 1.4, np.sin(r) * 1.4)
        else:
            new = (np.cos(r) * 1.4, np.sin(r) * 1.4)
        
        # plot circles
        circ = patches.Circle(new, radius=0.1, transform=ax.transData, edgecolor='k', facecolor=df[sequence[i]])
        circ.set_zorder(2)  # level in front of lines
        ax.add_patch(circ)
        
        # check if N- or C-terminus and add subscript, then plot AA letter
        if i == 0:
            ax.text(new[0], new[1], sequence[i] + '$_N$', va='center', ha='center', transform=ax.transData,
                    size=12, color=dt[sequence[i]], fontweight='bold')
        elif i == len(sequence) - 1:
            ax.text(new[0], new[1], sequence[i] + '$_C$', va='center', ha='center', transform=ax.transData,
                    size=12, color=dt[sequence[i]], fontweight='bold')
        else:
            ax.text(new[0], new[1], sequence[i], va='center', ha='center', transform=ax.transData,
                    size=12, color=dt[sequence[i]], fontweight='bold')
        
        eb = d_eisberg[sequence[i]][0]  # eisenberg value for this AA
        hm.append([eb * new[0], eb * new[1]])  # save eisenberg hydrophobicity vector value to later calculate HM
        
        old = (np.cos(r), np.sin(r))  # save as previous coordinates
    
    # draw hydrophobic moment arrow if moment option
    if moment:
        v_hm = np.sum(np.array(hm), 0)
        x = .0333 * v_hm[0]
        y = .0333 * v_hm[1]
        ax.arrow(0., 0., x, y, head_width=0.04, head_length=0.03, transform=ax.transData,
                 color='k', linewidth=6.)
        desc = calculate_hydrophobic_moment(sequence)
        if abs(x) < 0.2 and y > 0.:  # right positioning of HM text so arrow does not cover it
            z = -0.2
        else:
            z = 0.2
        plt.text(0., z, str(round(desc, 3)), fontdict={'fontsize': 20, 'fontweight': 'bold',
                                                                        'ha': 'center'})
    
    # plot shape
    if len(sequence) < 19:
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
    else:
        ax.set_xlim(-1.4, 1.4)
        ax.set_ylim(-1.4, 1.4)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    cur_axes = plt.gca()
    cur_axes.axes.get_xaxis().set_visible(False)
    cur_axes.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    
    if seq:
        plt.title(sequence, fontweight='bold', fontsize=20)
    
    # show or save plot
    if filename:
        plt.savefig(filename, dpi=150)
    else:
        pass



# Load Hemolysis model
with open('hemolysis_model.pkl', 'rb') as f:
    model_hemo = pickle.load(f)

 # Load Hemolysis vectorizer
with open('hemolysis_vectorizer.pkl', 'rb') as f:
    vectorizer_hemo = pickle.load(f)    

class DataFrameViewer(tk.Toplevel):
    def __init__(self, parent, df):
        tk.Toplevel.__init__(self, parent)
        self.title('Data Viewer')

        # Override the "X" button behavior
        self.protocol("WM_DELETE_WINDOW", self.disable_close)

        for column in df.select_dtypes(include=[np.number]).columns:
            df[column] = df[column].apply(lambda x: round_to_sig_figs(x, 2))
        self.images = []

        self.tree = tk.ttk.Treeview(self)
        self.tree["columns"] = list(df.columns)
        self.tree['show'] = 'headings'
        for column in self.tree["columns"]:
            self.tree.heading(column, text=column)
            self.tree.column(column, width=100)  # adjust width as needed
        for index, row in df.iterrows():
            self.tree.insert("", "end", values=list(row))

        self.tree.pack(fill='both', expand=True)

        # Add buttons at the bottom
        tk.Button(self, text="Mutagenesis", command=self.mutagenesis).pack(side=tk.LEFT)
        tk.Button(self, text="Optimization", command=self.optimization).pack(side=tk.LEFT)
        tk.Button(self, text="Close", command=lambda: [self.destroy(), self.master.ask_more_sequences()]).pack(side=tk.RIGHT)

    def disable_close(self):
        # This function does nothing, effectively disabling the "X" button
        pass

    def mutagenesis(self):
        
        # Check if an item is selected in the tree
        if not self.tree.selection():
            # No item is selected, show a warning popup
            messagebox.showwarning("Warning", "You must select a peptide before attempting mutagenesis or optimization")
            return

        # Get the selected item
        selected_item = self.tree.selection()[0]

        # Extract the peptide sequence from the selected item
        sequence = self.tree.item(selected_item)['values'][1]

        # If the sequence ends with 'X', remove it
        if sequence.endswith('X'):
            sequence = sequence[:-1]

        # Generate all possible point mutations
        mutated_sequences = []
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        for i in range(len(sequence)):
            for aa in amino_acids:
                mutated_sequences.append(sequence[:i] + aa + sequence[i+1:])

        # Load the AMPValidate model and vectorizer
        with open("AMPValidate.pkl", 'rb') as file:
            model_amp = pickle.load(file)
        with open("amp_validate_vectorizer.pkl", 'rb') as file:
            vectorizer_amp = pickle.load(file)

        # Load the Hemolysis model and vectorizer
        with open("hemolysis_model.pkl", 'rb') as file:
            model_hemo = pickle.load(file)
        with open("hemolysis_vectorizer.pkl", 'rb') as file:
            vectorizer_hemo = pickle.load(file)

        # Preprocess the mutated sequences and transform them
        sequences_split = [split_sequence(seq) for seq in mutated_sequences]
        X_amp = vectorizer_amp.transform(sequences_split)
        X_hemo = vectorizer_hemo.transform(sequences_split)

        # Make predictions
        probabilities_amp = model_amp.predict_proba(X_amp)[:, 1]  # take the second value from the output of predict_proba
        probabilities_hemo = model_hemo.predict_proba(X_hemo)[:, 1]  # take the second value from the output of predict_proba

        # Calculate the half-lives
        half_lives = [calc_half_life(seq) for seq in mutated_sequences]

        # Reshape the probabilities and half-lives into a 2D array for the heatmap
        prob_matrix_amp = np.array(probabilities_amp).reshape(len(sequence), len(amino_acids))
        prob_matrix_hemo = np.array(probabilities_hemo).reshape(len(sequence), len(amino_acids))
        half_life_matrix = np.array(half_lives).reshape(len(sequence), len(amino_acids))

        # Create a DataFrame for the heatmap
        heatmap_df_amp = pd.DataFrame(prob_matrix_amp, columns=list(amino_acids), index=list(sequence))
        heatmap_df_hemo = pd.DataFrame(prob_matrix_hemo, columns=list(amino_acids), index=list(sequence))
        heatmap_df_half_life = pd.DataFrame(half_life_matrix, columns=list(amino_acids), index=list(sequence))

        # Calculate the fitness scores
        w1, w2, w3 = 1, 1, 1
        fitness_scores = fitness(probabilities_amp, probabilities_hemo, half_lives, w1, w2, w3)

        # Reshape the fitness scores into a 2D array for the heatmap
        fitness_matrix = np.array(fitness_scores).reshape(len(sequence), len(amino_acids))

        # Create a DataFrame for the heatmap
        heatmap_df_fitness = pd.DataFrame(fitness_matrix, columns=list(amino_acids), index=list(sequence))

        # Create the heatmaps
        fig, axs = plt.subplots(2, 2, figsize=(12, 12))

        sns.heatmap(heatmap_df_amp, cmap='YlGnBu', ax=axs[0, 0])
        axs[0, 0].set_title('AMPValidate Probabilities for Point Mutations')

        sns.heatmap(heatmap_df_hemo, cmap='YlGnBu', ax=axs[0, 1])
        axs[0, 1].set_title('Hemolysis Probabilities for Point Mutations')

        sns.heatmap(heatmap_df_half_life, cmap='YlGnBu', ax=axs[1, 0])
        axs[1, 0].set_title('Half-lives for Point Mutations')

        sns.heatmap(heatmap_df_fitness, cmap='YlGnBu', ax=axs[1, 1])
        axs[1, 1].set_title('Fitness Scores for Point Mutations')

        plt.tight_layout()
        plt.show()

    def optimization(self):

        if not self.tree.selection():
        # No item is selected, show a warning popup
            messagebox.showwarning("Warning", "You must select a peptide before attempting mutagenesis or optimization")
            return
      
        # Get the selected item
        selected_item = self.tree.selection()[0]

        # Extract the peptide sequence from the selected item
        sequence = self.tree.item(selected_item)['values'][1]
        sequence = sequence.replace('C', 'S')

        # If the sequence ends with 'X', remove it
        if sequence.endswith('X'):
            sequence = sequence[:-1]

        # Create a new top-level window
        wait_window = tk.Toplevel(self)
        center_window(wait_window)
        wait_window.overrideredirect(True)

        # Set the size of the window
        window_width = 600  # Width in pixels
        window_height = 300  # Height in pixels
        wait_window.geometry(f"{window_width}x{window_height}")

        # Add a label with the main message
        label_text_main = "Please wait while I am optimizing your sequence."
        tk.Label(wait_window, text=label_text_main, wraplength=window_width - 20, font=("Helvetica", 14)).pack()
        
        # Add some vertical space
        tk.Label(wait_window, text="").pack(pady=10) # Adjust the pady value to change the space

        # Add a label with the ATTENTION message, in a larger and bold font
        label_text_attention = "ATTENTION: This optimizer only considers linear peptides without disulfide bonds, so all cysteine residues will be replaced by serine"
        tk.Label(wait_window, text=label_text_attention, wraplength=window_width - 20, font=("Helvetica", 18, "bold")).pack()


        # Add the logo
        logo = tk.PhotoImage(file="progress.png")
        logo_label = tk.Label(wait_window, image=logo)
        logo_label.pack()

        # Add the progress bar
        progress = ttk.Progressbar(wait_window, mode='indeterminate')
        progress.pack()
        progress.start()

        # Center the window
        wait_window.update()  # Update to get the correct window dimensions
        window_width = wait_window.winfo_width()
        window_height = wait_window.winfo_height()
        position_right = int(wait_window.winfo_screenwidth()/2 - window_width/2)
        position_down = int(wait_window.winfo_screenheight()/2 - window_height/2)
        wait_window.geometry("+{}+{}".format(position_right, position_down))


        # Initialize the population with the selected sequence and its mutated versions
        population = [sequence] + [mutate_sequence(sequence) for _ in range(POPULATION_SIZE - 1)]  # mutate_sequence is a function that applies random mutations to a sequence

        # Initialize the best fitness and best sequence
        best_fitness = -1
        best_sequence = None

        # Counter for generations with no improvement
        no_improvement_counter = 0

        for generation in range(NUM_GENERATIONS):
            # Evaluation
            fitness_scores = []
            for seq in population:
                prob_amp = predict_proba_amp(seq)  # You need to define this function.
                prob_hemo = predict_proba_hemo(seq)  # You need to define this function.
                half_life = calc_half_life(seq)
                fitness_score = fitness(prob_amp, prob_hemo, half_life, w1, w2, w3)
                fitness_scores.append(fitness_score)

            # Check for improvement
            max_fitness = max(fitness_scores)
            if max_fitness > best_fitness:
                best_fitness = max_fitness
                best_sequence = population[fitness_scores.index(max_fitness)]
                no_improvement_counter = 0  # reset the counter
            else:
                no_improvement_counter += 1

            # Check the stopping criterion
            if no_improvement_counter >= MAX_NO_IMPROVEMENT:
                break

            # Selection
            selected_individuals = random.choices(population, weights=fitness_scores, k=POPULATION_SIZE)

            # Crossover
            offspring = []
            for i in range(0, POPULATION_SIZE, 2):
                parent1 = selected_individuals[i]
                parent2 = selected_individuals[i+1]
                if random.random() < CROSSOVER_RATE:
                    crossover_point = random.randint(1, len(parent1) - 1)
                    child1 = parent1[:crossover_point] + parent2[crossover_point:]
                    child2 = parent2[:crossover_point] + parent1[crossover_point:]
                else:
                    child1, child2 = parent1, parent2
                offspring.append(child1)
                offspring.append(child2)

            # Mutation
            for i in range(POPULATION_SIZE):
                if random.random() < MUTATION_RATE:
                    mutation_point = random.randint(0, len(offspring[i]) - 1)
                    # Exclude cysteine ('C') from the possible amino acids
                    new_amino_acid = random.choice([aa for aa in amino_acids if aa != 'C'])
                    offspring[i] = offspring[i][:mutation_point] + new_amino_acid + offspring[i][mutation_point+1:]

            # Replacement
            population = offspring

        # Print the best sequence and its fitness score
        print(f"Best sequence: {best_sequence} with fitness score: {best_fitness}")

        # Close the waiting popup
        wait_window.destroy()

        # Generate and save the helical wheel plots without displaying them
        helical_wheel(sequence, moment=True, filename='original_helix.png')
        helical_wheel(best_sequence, moment=True, filename='optimized_helix.png')

        # Display the results
        results_window = tk.Toplevel(self.master)
        results_window.overrideredirect(True)

        # Add the original and optimized sequences (and their properties)
        prob_amp_orig = predict_proba_amp(sequence)
        prob_hemo_orig = predict_proba_hemo(sequence)
        half_life_orig = calc_half_life(sequence)
        prob_amp_opt = predict_proba_amp(best_sequence)
        prob_hemo_opt = predict_proba_hemo(best_sequence)
        half_life_opt = calc_half_life(best_sequence)
        original_sequence_label = tk.Label(results_window, text=f'Original sequence: {sequence} \n Probability of antimicrobial: {prob_amp_orig} \n Probability of hemolytic: {prob_hemo_orig} \n Half-life: {half_life_orig}')
        original_sequence_label.pack()
        optimized_sequence_label = tk.Label(results_window, text=f'Optimized sequence: {best_sequence} \n Probability of antimicrobial: {prob_amp_opt} \n Probability of hemolytic: {prob_hemo_opt} \n Half-life: {half_life_opt}')
        optimized_sequence_label.pack()

        # Load the helix wheel images
        original_helix_img = Image.open("original_helix.png")
        optimized_helix_img = Image.open("optimized_helix.png")

        # Convert images to Tkinter-compatible format
        tk_img_orig = ImageTk.PhotoImage(original_helix_img)
        tk_img_opt = ImageTk.PhotoImage(optimized_helix_img)

        # Create labels and add the images
        label_orig = tk.Label(results_window, image=tk_img_orig)
        label_orig.pack(side=tk.LEFT)
        label_opt = tk.Label(results_window, image=tk_img_opt)
        label_opt.pack(side=tk.RIGHT)

        # Close the images
        original_helix_img.close()
        optimized_helix_img.close()

        # Close the results window when the user clicks the 'Close' button
        close_button = tk.Button(results_window, text="Close", command=results_window.destroy)
        close_button.pack()

        results_window.mainloop()


class Application(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.title("PyAMPA")
        self.geometry("1200x750")

        # Center the main window
        center_window(self)

        # Load and resize the logo
        image = Image.open("logo.gif")  # replace 'logo.gif' with the path to your logo file
        image = image.resize((500, 200), Image.LANCZOS)  # specify the size in pixels
        self.logo = ImageTk.PhotoImage(image)

        tk.Label(self, image=self.logo).pack()

        # Welcome message
        tk.Label(self, text="Welcome to PyAmpa! Please upload your sequences in fasta file", font=("Helvetica", 20, 'bold')).pack()

        # Add a button to select the input file
        self.input_button = tk.Button(self, text="Input", command=self.select_input_file)
        self.input_button.pack()

        # Add a label to display the path of the selected file
        self.file_path_label = tk.Label(self, text="")
        self.file_path_label.pack(pady=20)  # Increase vertical padding to add distance

    def select_input_file(self):
        self.input_file = filedialog.askopenfilename(title="Select fasta file",
                                                     filetypes=(("fasta files", "*.fasta"), ("all files", "*.*")))
        self.file_path_label.config(text=self.input_file)

        self.output_message_label = tk.Label(self, text="Please select the output directory. All files will be saved in this directory", font=("Helvetica", 20, 'bold'))
        self.output_message_label.pack(pady=20)  # Increase vertical padding to add distance

        # Add a button to select the output directory
        self.output_button = tk.Button(self, text="Output", command=self.select_output_directory)
        self.output_button.pack()

        # Add a label to display the selected output directory
        self.output_directory_label = tk.Label(self, text="")
        self.output_directory_label.pack()

    def select_output_directory(self):
        self.output_directory = filedialog.askdirectory(title="Select output directory")
        self.output_directory_label.config(text=self.output_directory)

        # Wait one second before launching the popup window
        self.after(10, self.launch_popup)

    def launch_popup(self):
        # Create a popup window for running pyampa.py
        self.popup = tk.Toplevel(self)
        self.popup.overrideredirect(True)
        self.popup.geometry("600x300")  # Adjust the dimensions
        self.popup.geometry(f"+{self.winfo_x() + (self.winfo_width() - 600) // 2}+{self.winfo_y() + (self.winfo_height() - 300) // 2}")

        # Load and resize the man running logo
        man_running_image = Image.open("man_running.png")  # replace 'man_running.png' with the path to your image file
        man_running_image = man_running_image.resize((100, 100), Image.LANCZOS)  # specify the size in pixels
        self.man_running_logo = ImageTk.PhotoImage(man_running_image)

        # Add a label to display the logo
        tk.Label(self.popup, image=self.man_running_logo).pack()

        # Add a label to display a message
        self.run_label = tk.Label(self.popup, text="Running the prediction, please wait until completion", font=("Helvetica", 14))
        self.run_label.pack(pady=10)

        # Add a progress bar
        self.progress = Progressbar(self.popup, orient=tk.HORIZONTAL, length=300, mode='indeterminate')
        self.progress.pack(pady=20)
        self.progress.start()

        # Run pyampa.py
        self.run_pyampa()

    def run_pyampa(self):
        # Run pyampa.py
        subprocess.run(["python", "pyampa.py", "-i", self.input_file, "-o", self.output_directory], check=True)

        # Wait one second before closing the popup window
        self.after(1000, self.close_popup)

    def close_popup(self):
        # Stop the progress bar
        self.progress.stop()

        # Close the popup window
        self.popup.destroy()

        # Display results summary
        answer = messagebox.askquestion("Summary", "Your prediction is completed. Do you want to print a summary of your results?")
        if answer == "yes":
            self.create_summary()
        else:
            # If the user does not want to print a summary, ask if they want to run AMPValidate
            self.ask_amp_validate()

    def create_summary(self):
        # Read results.csv
        results_file = f"{self.output_directory}/results.csv"
        df = pd.read_csv(results_file)

        # Create a single figure with multiple subplots
        plt.figure(figsize=(10, 8))

        # Histogram of peptide lengths
        plt.subplot(2, 2, 1)
        df['Peptide Length'] = df['Peptide'].apply(len)
        plt.hist(df['Peptide Length'], bins=20, color='skyblue')
        plt.xlabel('Peptide Length')
        plt.ylabel('Frequency')
        plt.title('Histogram of Peptide Lengths')

        # Amino acid composition diagram
        plt.subplot(2, 2, 2)
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        composition = df['Peptide'].apply(lambda x: [x.count(aa) for aa in amino_acids])
        composition = pd.DataFrame(composition.tolist(), columns=list(amino_acids))
        composition_sum = composition.sum()
        composition_sum.plot(kind='bar', color='lightcoral')
        plt.xlabel('Amino Acid')
        plt.ylabel('Frequency')
        plt.title('Amino Acid Composition of Predicted Peptides')

        # Histogram of antimicrobial indices
        plt.subplot(2, 2, 3)
        plt.hist(df['Average Torrent Index'], bins=20, color='lightgreen')
        plt.xlabel('Antimicrobial Index')
        plt.ylabel('Frequency')
        plt.title('Histogram of Antimicrobial Indices')

        # Remove empty plot for better layout
        # Histogram of hydrophobic moments
        plt.subplot(2, 2, 4)
        df['Hydrophobic Moment'] = df['Peptide'].apply(calculate_hydrophobic_moment)
        plt.hist(df['Hydrophobic Moment'], bins=20, color='purple')
        plt.xlabel('Hydrophobic Moment')
        plt.ylabel('Frequency')
        plt.title('Histogram of Hydrophobic Moments')

        # Save the plots in the output directory
        plt.tight_layout()
        plt.savefig(f"{self.output_directory}/summary_plots.png")

        # Show the plots
        plt.show()

        # After the user closes the plots window, ask if they want to run AMPValidate
        self.ask_amp_validate()

    def ask_amp_validate(self):
        # Ask if the user wants to run AMPValidate
        answer = messagebox.askquestion("AMPValidate", "Do you want to run AMPValidate on your predicted sequences?")
        if answer == "yes":
            self.run_amp_validate()
        else:
            self.ask_more_sequences()

    def show_dataframe(self, df):
        DataFrameViewer(self, df)

    def run_amp_validate(self):
        # Load the model and the vectorizer
        with open("AMPValidate.pkl", 'rb') as file:
            model = pickle.load(file)
        with open("amp_validate_vectorizer.pkl", 'rb') as file:
            vectorizer = pickle.load(file)

        # Read results.csv
        results_file = f"{self.output_directory}/results.csv"
        df = pd.read_csv(results_file)

        # Preprocess the peptide sequences and transform them
        sequences_split = df['Peptide'].apply(split_sequence)
        X = vectorizer.transform(sequences_split)

        # Make predictions
        df['AMPValidate Probability'] = model.predict_proba(X)[:, 1]  # take the second value from the output of predict_proba

        # Save the predictions in the results.csv file
        df.to_csv(results_file, index=False)



        # Load the results.csv file
        results_df = pd.read_csv(results_file)

        # Load the trained model
        with open('activities_model.pkl', 'rb') as f:
            model = pickle.load(f)

        # Load the fitted CountVectorizer
        with open('activities_vectorizer.pkl', 'rb') as f:
            vectorizer = pickle.load(f)

        # Load the LabelEncoder
        with open('label_encoder.pkl', 'rb') as f:
            label_encoder = pickle.load(f)

        # Add an 'X' at the end of each peptide sequence to amidate at the C.terminus
        results_df['Peptide'] = results_df['Peptide'].apply(lambda x: x + 'X')

        # Preprocess the peptide sequences
        results_df['sequence_split'] = results_df['Peptide'].apply(split_sequence)

        # Count the occurrences of each subsequence using the fitted CountVectorizer
        X_results = vectorizer.transform(results_df['sequence_split'])

        # Create a new DataFrame with the counts
        counts_results_df = pd.DataFrame(X_results.toarray(), columns=vectorizer.get_feature_names_out())

        # Make predictions for each bacterium
        predictions = []
        for bacterium in range(len(label_encoder.classes_)):
            counts_results_df['bacterium'] = bacterium
            prediction = model.predict(counts_results_df)
            predictions.append(prediction)

        # Transpose the list of predictions and convert it to a DataFrame
        predictions_df = pd.DataFrame(list(map(list, zip(*predictions))), columns=label_encoder.classes_)

        # Add the predictions to the results DataFrame
        results_df = pd.concat([results_df, predictions_df], axis=1)

        # Drop the 'sequence_split' column
        results_df = results_df.drop(columns=['sequence_split'])

        # Undo the log transformation for the bacteria predictions
        for bacterium in label_encoder.classes_:
            results_df[bacterium] = np.exp(results_df[bacterium])

        # Remove the `X` residue 
        results_df['Peptide'] = results_df['Peptide'].str.replace('X', '')

        # Save the updated DataFrame to the results.csv file
        results_df.to_csv(results_file, index=False)




        # Load the model and the vectorizer for hemolysis
        with open("hemolysis_model.pkl", 'rb') as file:
            hemo_model = pickle.load(file)
        with open("hemolysis_vectorizer.pkl", 'rb') as file:
            hemo_vectorizer = pickle.load(file)

        # Load the model and the vectorizer for cell-penetrating capacity
        with open("cpp_model.pkl", 'rb') as file:
            cpp_model = pickle.load(file)
        with open("cpp_vectorizer.pkl", 'rb') as file:
            cpp_vectorizer = pickle.load(file)

        # Load the model and the vectorizer for toxicity
        with open("tox_model.pkl", 'rb') as file:
            tox_model = pickle.load(file)
        with open("tox_vectorizer.pkl", 'rb') as file:
            tox_vectorizer = pickle.load(file)

        # Read results.csv
        results_file = f"{self.output_directory}/results.csv"
        results_df = pd.read_csv(results_file)

        # Preprocess the peptide sequences and transform them
        sequences_split = results_df['Peptide'].apply(split_sequence)

        # Make predictions for hemolysis
        X_hemo = hemo_vectorizer.transform(sequences_split)
        results_df['Hemolytic Probability'] = hemo_model.predict_proba(X_hemo)[:, 1]  # take the second value from the output of predict_proba

        # Make predictions for cell-penetrating capacity
        X_cpp = cpp_vectorizer.transform(sequences_split)
        results_df['Cell-penetrating Probability'] = cpp_model.predict_proba(X_cpp)[:, 1]  # take the second value from the output of predict_proba

        # Make predictions for toxicity
        X_tox = tox_vectorizer.transform(sequences_split)
        results_df['Toxic Probability'] = tox_model.predict_proba(X_tox)[:, 1]  # take the second value from the output of predict_proba

        # Save the predictions in the results.csv file
        results_df.to_csv(results_file, index=False)


        # Select only the columns to display
        results_df = results_df[['Name', 'Peptide', 'Average Torrent Index', 'AMPValidate Probability', 'Hemolytic Probability', 'Toxic Probability', 'Cell-penetrating Probability']]

        # Rename the columns (if necessary)
        results_df.columns = ['Name', 'Peptide', 'Avg Torrent Index', 'AMP Validate Prob', 'Hemolytic Prob', 'Toxic Prob', 'Cell-penetrating Prob']



        # Load and resize the validation logo
        validate_image = Image.open("validate.png")  # replace 'validate.png' with the path to your image file
        validate_image = validate_image.resize((100, 100), Image.LANCZOS)  # specify the size in pixels
        self.validate_logo = ImageTk.PhotoImage(validate_image)

        # Create a new window to display the logo and the message
        validate_window = tk.Toplevel(self)
        validate_window.overrideredirect(True)
        tk.Label(validate_window, image=self.validate_logo).pack()
        tk.Label(validate_window, text="Your sequences have been validated. Happy testing!").pack()
        
        # Call the show_dataframe function when the OK button is clicked
        ok_button_command = lambda: [validate_window.destroy(), self.show_dataframe(results_df)]
        tk.Button(validate_window, text="OK", command=ok_button_command).pack()


        # Center the window on the screen
        validate_window.update()  # update window size before centering
        width = validate_window.winfo_width()
        height = validate_window.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2) -150
        validate_window.geometry(f'{width}x{height}+{x}+{y}')

        # Wait until the user clicks "OK"
        self.wait_window(validate_window)

    def ask_more_sequences(self):
        # Load and resize the continue logo
        continue_image = Image.open("continue.png")  # replace 'continue.png' with the path to your image file
        continue_image = continue_image.resize((100, 100), Image.LANCZOS)  # specify the size in pixels
        self.continue_logo = ImageTk.PhotoImage(continue_image)

        # Create a new window to display the logo and ask if the user wants to process more sequences
        continue_window = tk.Toplevel(self)
        continue_window.title("More sequences?")
        tk.Label(continue_window, image=self.continue_logo).pack()
        tk.Label(continue_window, text="Happy testing! Do you want to process more sequences?").pack()
        yes_button = tk.Button(continue_window, text="Yes", command=lambda: [continue_window.destroy(), self.restart()])
        yes_button.pack(side=tk.LEFT)
        no_button = tk.Button(continue_window, text="No", command=lambda: [continue_window.destroy(), self.show_citation()])
        no_button.pack(side=tk.RIGHT)

        # Center the window on the screen
        continue_window.update()  # update window size before centering
        width = continue_window.winfo_width()
        height = continue_window.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        continue_window.geometry(f'{width}x{height}+{x}+{y}')

    def show_citation(self):
        # Load and resize the citation logo
        citation_image = Image.open("citation.png")  # replace 'citation.png' with the path to your image file
        citation_image = citation_image.resize((100, 100), Image.LANCZOS)  # specify the size in pixels
        self.citation_logo = ImageTk.PhotoImage(citation_image)

        # Create a new window to display the logo and show the citation message
        citation_window = tk.Toplevel(self)
        citation_window.title("Citation")
        tk.Label(citation_window, image=self.citation_logo).pack()
        tk.Label(citation_window, text="Thanks for using pyAMPA. If you find this program useful, please cite us.").pack()
        finish_button = tk.Button(citation_window, text="Finish", command=lambda: os._exit(0))
        finish_button.pack()

        # Center the window on the screen
        citation_window.update()  # update window size before centering
        width = citation_window.winfo_width()
        height = citation_window.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        citation_window.geometry(f'{width}x{height}+{x}+{y}')


    def restart(self):
        # Restart the script from the beginning
        self.input_button.destroy()
        self.file_path_label.destroy()
        self.output_message_label.destroy()
        self.output_button.destroy()
        self.output_directory_label.destroy()

         # Reset the file paths
        self.input_file_path = None
        self.output_directory_path = None
        
        # Add a button to select the input file
        self.input_button = tk.Button(self, text="Input", command=self.select_input_file)
        self.input_button.pack()

        # Add a label to display the path of the selected file
        self.file_path_label = tk.Label(self, text="")
        self.file_path_label.pack(pady=20)  # Increase vertical padding to add distance


# Start the application
app = Application()
app.mainloop()
