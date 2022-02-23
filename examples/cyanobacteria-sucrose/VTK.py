import vtk
colors = vtk.vtkNamedColors()
atom = './Run_1_1/atom_0.vtu'
box = './Run_1_1/atom_0_boundingBox.vtr'
    # Read the source file.
    # Create the reader for the data.
# print("Loading ", filename)
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(atom)
reader.Update()
output = reader.GetOutput()
diameters = output.GetPointData().GetArray('diameter')
types = output.GetPointData().GetArray('type')
points = output.GetPointData()


# Glyph the points
sphere = vtk.vtkSphereSource()

# sphere.SetRadius(diameters)
# # Set up the glyph filter
glyph = vtk.vtkGlyph3D()
# glyph.SetInputConnection(output.GetOutputPort())
# glyph.SetSourceConnection(sphere.GetOutputPort())


# # Create the mapper that creates graphics elements
mapper = vtk.vtkDataSetMapper()
mapper.SetInputConnection(reader.GetOutputPort())

# # Create the Actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)
# # show the edges of the image grid
actor.GetProperty().SetRepresentationToWireframe()
# actor.GetProperty().SetColor(colors.GetColor3d("DarkSalmon"))

# # Create the Renderer
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.ResetCamera()
# renderer.SetBackground(colors.GetColor3d("Silver"))

# # Create the RendererWindow
renderer_window = vtk.vtkRenderWindow()
renderer_window.AddRenderer(renderer)

# # Create the RendererWindowInteractor and display the vti file
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderer_window)
interactor.Initialize()
interactor.Start()
