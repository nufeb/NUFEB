import vtk
colors = vtk.vtkNamedColors()
filename = './Sucrose_0.0/atom_0.vtu'
    # Read the source file.
    # Create the reader for the data.
print("Loading ", filename)
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(filename)
reader.Update()



# Glyph the points
sphere = vtk.vtkSphereSource()
# sphere.SetPhiResolution(21)
# sphere.SetThetaResolution(21)
sphere.SetRadius(.08)



# # Create the mapper that creates graphics elements
# mapper = vtk.vtkDataSetMapper()
# mapper.SetInputConnection(reader.GetOutputPort())

# # Create the Actor
# actor = vtk.vtkActor()
# actor.SetMapper(mapper)
# # show the edges of the image grid
# actor.GetProperty().SetRepresentationToWireframe()
# actor.GetProperty().SetColor(colors.GetColor3d("DarkSalmon"))

# # Create the Renderer
# renderer = vtk.vtkRenderer()
# renderer.AddActor(actor)
# renderer.ResetCamera()
# renderer.SetBackground(colors.GetColor3d("Silver"))

# # Create the RendererWindow
# renderer_window = vtk.vtkRenderWindow()
# renderer_window.AddRenderer(renderer)

# # Create the RendererWindowInteractor and display the vti file
# interactor = vtk.vtkRenderWindowInteractor()
# interactor.SetRenderWindow(renderer_window)
# interactor.Initialize()
# interactor.Start()