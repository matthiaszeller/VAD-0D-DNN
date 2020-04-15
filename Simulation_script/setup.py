
# =================== CONFIGURATION FILE

# ============ DATASET GENERATION

# Number of samples for dataset generation
numberofsamples = 10000

# Original modelica file
file_path="/media/maousi/Data/Documents/Programmation/git/vad-0d-dnn/modelica/original/Mathcard.mo"

# Output folder path for the raw dataset
output_folder="/media/maousi/Data/tmp/simulation_noLVAD"

# SIMULATION WITH/WITHOUT LVAD
# Warning: this applies for dataset generation and for DNN testing
SIMULATION_LVAD = False

# Verbose and debug mode
DEBUG_MODE = True
def q(msg):
    if DEBUG_MODE:
        print("<Debug msg> " + msg)


# ============ DNN TESTING

# Location of the DNN training files
dnn_folder = "/media/maousi/Data/tmp/simulations_2020_03_21/dnn_16neurons"

# Output folder path for the DNN testing
output_folder_DNN_test = "/media/maousi/Data/tmp/simulations_2020_03_21/dnn_16neu_test"
