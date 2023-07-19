DEFAULT_CONFIG_FILE_PATH = 'configuration_help.json'
DEFAULT_BIN_FILE_PATH = '/home/ruben/PycharmProjects/DensityPy/bin'
DEFAULT_CONFIG_CONTENT = {
    "Project settings": {
        "project_name": "NMA",
        "XYZ Molecule Geometry": "NMA.xyz",
        # "Basis": "ANO-L-MB",
        "Number of States": 4,
        "List_of_Active_Orbitals": [ 18, 19, 20, 21],
        "Molcas Output Directory": "NMA_output",
    },
    "GRID settings": {
        "Number of Points X axis": 15,
        "Number of Points Y axis": 15,
        "Number of Points Z axis": 15,
        "X MIN": -5,
        "X MAX": 6,
        "Y MIN": -4,
        "Y MAX": 5,
        "Z MIN": -4,
        "Z MAX": 4,
        "Step Size": 0.1,
        "Boundary": 2
    },
    "Charge Migration settings": {
        "Experiment Directory": "sim",
        "Field File": "field_file",
        "Number of Times": 1001,
        "Min Time": -4000,
        "Max Time": 5000,
        "Bath Temperature": 3273.75,
        "Dephasing Factor": 0.001,
        "Relaxation Factor": 0.001
    },
    "Pump Pulses settings": {
        "Type of Pulse": "G",
        "start_time": 0.000,
        "Pump Central Frequency": 0.4537,
        "Pump Periods": 5,
        "Pump Phase": 0,
        "Pump Intensity": 0.00001,
        "Pump Polarization": [90, 90]
    },
    "Probe Pulses settings": {
        "Type of Pulse": "G",
        "Time Delay Start": -2000,
        "Time Delay Stop": 3000,
        "Number Of PP": 101,
        "Probe Central Frequency": 0.07,
        "Probe Periods": 5,
        "Probe Phase": 0,
        "Probe Intensity": 0.00001,
        "Probe Polarization": [90, 90]
    },
    "Charge Migration FT settings": {
        "Number of Omegas": 101,
        "Min Omega": 0.2,
        "Max Omega": 0.65,
        "Number of TauOmegas": 101,
        "Min TauOmega": -0.20,
        "Max TauOmega": 0.20,
        "FT TimeStep": 3200,
        "FT WidthStep": 400
    }
}
