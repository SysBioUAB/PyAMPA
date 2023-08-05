
# PyAMPA: A new tool for Antimicrobial Peptide Discovery and Optimization

## Introduction

PyAMPA is a complete new implementation of the AMPA algorithm (https://tcoffee.crg.eu/apps/ampa/do), originally designed to predict antimicrobial sequences in proteins. With PyAMPA you can now predict antimicrobial sequences from entire proteomes and has new features such as optimization, mutagenesis and interactive tables for visualizing data.

**If you found this work useful, please cite us:**

**--PyAMPA-- Coming soon.**

**--Original AMPA-- Marc Torrent, Paolo Di Tommaso, David Pulido, M. Victòria Nogués, Cedric Notredame, Ester Boix, David Andreu: AMPA: an automated web server for prediction of protein antimicrobial regions (2011) Bioinformatics 28(1):130-131.**


## Installation

Follow the instructions below to set up and run the application.

### Prerequisites

- Python 3 (tested with Python 3.11)
  Alreday installed by default on macOS and Linux OS but if using Windows you may need to install Python 3
  Please check your current version of Python

  ```bash
  python --version
  ```

  If your version is < 3.11, updating Python is suggested

- Tkinter for the GUI

  To install Tkinter in Ubuntu Linux:

  ```bash
  sudo apt install python3-tk
  ```
  
  To install Tkinter:

  ```bash
  pip install tk
  ```
  
- Additional libraries: NumPy, Pandas, Matplotlib, Seaborn, Scikit-learn, PIL (Pillow). You can install the required dependencies using the following command:

```bash
pip install numpy pandas matplotlib seaborn scikit-learn pillow
```

### Downloading the Application

Download the `main.py` file and all other files to your local machine.

## Running the Application

After installing the required dependencies and downloading the application, follow the steps below to run it:

1. Navigate to the directory containing the `main.py` file.
2. Run the following command:

```bash
python main.py
```

The application should open, and you can interact with it through the graphical user interface.

## Usage

Within the application, you can:

- Predict antimicrobial peptides from entire proteomes.
- Select peptide sequences for optimization or mutagenesis.
- Perform various analyses on peptide sequences and visualize the results using heatmaps and other graphical representations.

## Troubleshooting

If you encounter any issues or have questions about specific functionalities, please refer to the code documentation or contact the development team.

## License

Include information about the license if applicable.

## Contact

Provide contact details if users have further questions or need support.

---
