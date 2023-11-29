import csv
import vtk


def read_csv(file_name, fields):
    data_points = []

    with open(file_name, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            print(row)
            exit()
            data = tuple(row[field] for field in fields)
            data_points.append(data)

    print(data_points)
    exit()
    return data_points


def write_vtk(data_points, fields, output_file):
    points = vtk.vtkPoints()

    field_arrays = {field: vtk.vtkDoubleArray() for field in fields[3:]}
    for field in field_arrays:
        field_arrays[field].SetName(field)

    for data in data_points:
        x, y, z = map(float, data[:3])
        points.InsertNextPoint(x, y, z)
        for i, field in enumerate(fields[3:]):
            field_arrays[field].InsertNextValue(float(data[i + 3]))

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    for field in field_arrays:
        polydata.GetPointData().AddArray(field_arrays[field])

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(polydata)
    writer.Write()


if __name__ == "__main__":
    csv_file = "/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/sim/ChargeDensity/ChDenSimPP-4000.0/ChDen1.4014.csv"
    vtk_file = "output.vtk"

    # Define default fields
    default_fields = ["x", "y", "z", "ChargeDensity"]

    # Optionally extend this list based on user needs:
    # Example: ["x", "y", "z", "ChargeDensity", "Atom_O_ChargeDensity", ...]

    data_points = read_csv(csv_file, default_fields)
    write_vtk(data_points, default_fields, vtk_file)
