DEFAULT_BIN_FILE_PATH = '/home/ruben/PycharmProjects/DensityPy/bin'
DEFAULT_CONFIG_CONTENT = {
    "Project settings": {
        "project_name": "NMA",
        "XYZ Molecule Geometry": "NMA.xyz",
        "Number of States": 4,
        "List_of_Active_Orbitals": [
            18,
            19,
            20,
            21
        ],
        "Molcas Output Directory": "NMA_output",
        "Experiment Directory": "example_sim"
    },
    "GRID settings": {
        "Number of Points X axis": 100,
        "Number of Points Y axis": 100,
        "Number of Points Z axis": 100,
        "X MIN": -6,
        "X MAX": 7,
        "Y MIN": -5,
        "Y MAX": 6,
        "Z MIN": -5,
        "Z MAX": 5,
        "Step Size": 0.008,
        "Boundary": 6
    },
    "Charge Migration settings": {
        "Field File": "field_file",
        "Number of Times": 9001,
        "Min Time": -5000,
        "Max Time": 13000,
        "Bath Temperature": 3000,
        "Dephasing Factor": 0.0012,
        "Relaxation Factor": 0.0048
    },
    "Pump Pulses settings": {
        "Type of Pulse": "G",
        "start_time": 0,
        "Pump Central Frequency": 0.410,
        "Pump Periods": 2,
        "Pump Phase": 0,
        "Pump Intensity": 1.2e-1,
        "Pump Polarization": [
            90,
            90
        ]
    },
    "Probe Pulses settings": {
        "Type of Pulse": "G",
        "Time Delay Start": -200,
        "Time Delay Stop": 800,
        "Number Of PP": 120,
        "Time Delay Weight Factor": 0.5,
        "Probe Central Frequency": 0.073,
        "Probe Periods": 2,
        "Probe Phase": 0,
        "Probe Intensity": 1.2e-01,
        "Probe Polarization": [
            90,
            90
        ]
    },
    "Charge Migration FT settings": {
        "Number of Omegas": 120,
        "Min Omega": 0.18,
        "Max Omega": 0.68,
        "Number of TauOmegas": 120,
        "Min TauOmega": -0.25,
        "Max TauOmega": 0.25,
        "FT TimeStep": 5000,
        "FT WidthStep": 350
    }
}
