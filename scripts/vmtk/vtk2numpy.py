from paraview.simple import *


reader = OpenDataFile("tmp/Cross_Sectional_Area_vtp_file.vtp")
UpdatePipeline()
dataInfo = reader.GetDataInformation()

# print(reader.CellData['CenterlineSectionArea'].GetNumberOfComponents())
SaveData("Initialarea.csv", reader, Precision=5, FieldAssociation='Cell Data')
