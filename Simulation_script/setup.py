

# SETUP FILE
# Configure file paths

# SIMULATION WITH/WITHOUT LVAD
SIMULATION_LVAD = True

# Verbose and debug mode
DEBUG_MODE = True
def q(msg):
    if DEBUG_MODE:
        print("<Debug msg> " + msg)

# Original modelica file
filepath="/media/maousi/Data/Documents/Programmation/git/vad-0d-dnn/modelica/original/"

# Output folder path
outputfolder="/media/maousi/Data/tmp/OM_Simulation"
