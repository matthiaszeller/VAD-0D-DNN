

# SETUP FILE
# Configure file paths

numberofsamples = 10000

# SIMULATION WITH/WITHOUT LVAD
SIMULATION_LVAD = True

# Verbose and debug mode
DEBUG_MODE = True
def q(msg):
    if DEBUG_MODE:
        print("<Debug msg> " + msg)

# Original modelica file
file_path="/media/maousi/Data/Documents/Programmation/git/vad-0d-dnn/modelica/original/Mathcard.mo"

# Output folder path
output_folder="/media/maousi/Data/tmp/OM_SimulationNew"
