DEFAULT_CONFIG_FILE_PATH = 'configuration_help.json'
DEFAULT_CONFIG_CONTENT = {
    "Project settings": {
        "Project_Name": "NMA",
        "XYZ Molecule Geometry": "NMA.xyz",
        "Basis": "ANO-L-MB",
        "Number of States": "5",
        "List_of_Active_Orbitals": "17 18 19 20 21",
        "Molcas Output Directory": "NMA_output"
    },
    "GRID settings": {
        "Number of Points, X axis": "15",
        "Number of Points, Y axis": "15",
        "Number of Points, Z axis": "15",
        "X MIN": "-5",
        "X MAX": "6",
        "Y MIN": "-4",
        "Y MAX": "5",
        "Z MIN": "-4",
        "Z MAX": "4",
        "Step Size": "0.1",
        "Boundary": "2"
    },
    "Charge Migration settings": {
        "Output Directory": "sim",
        "Field File": "Field_File",
        "Number of Times": "1001",
        "Min Time": "-4000",
        "Max Time": "5000",
        "Bath Temperature": "3273.75",
        "Dephasing Factor": "0.001",
        "Relaxation Factor": "0.001"
    },
    "Pulses settings (Pump)": {
        "Type of Pulse": "G",
        "Start_Time": "0.000",
        "Pump Central Frequency": "0.4537",
        "Pump Periods": "5",
        "Pump Phase": "0",
        "Pump Intensity": "0.00001",
        "Pump Polarization": "90 90"
    },
    "Pulses settings (Probe)": {
        "Type of Pulse": "G",
        "Time Delay Start": "-2000",
        "Time Delay Stop": "3000",
        "Number Of PP": "101",
        "Probe Central Frequency": "0.07",
        "Probe Periods": "5",
        "Probe Phase": "0",
        "Probe Intensity": "0.00001",
        "Probe Polarization": "90 90"
    },
    "Charge Migration FT settings": {
        "Number of Omegas": "101",
        "Min Omega": "0.2",
        "Max Omega": "0.65",
        "Number of TauOmegas": "101",
        "Min TauOmega": "-0.20",
        "Max TauOmega": "0.20",
        "TimeStep (FT)": "3200",
        "WidthStep (FT)": "400"
    }
}
