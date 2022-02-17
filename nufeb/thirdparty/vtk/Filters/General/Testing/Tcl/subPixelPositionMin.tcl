package require vtk

# Test sub pixel positioning (A round about way of getting an iso surface.)
# See cubed sphere for the surface before sub pixel poisitioning.

vtkSphere sphere
    sphere SetCenter 1 1 1
    sphere SetRadius 0.9

vtkSampleFunction sample
    sample SetImplicitFunction sphere
    sample SetModelBounds  0 2 0 2 0 2
    sample SetSampleDimensions 30 30 30
    sample ComputeNormalsOff
    sample Update


vtkThreshold threshold1
    threshold1 SetInputConnection [sample GetOutputPort]
    threshold1 ThresholdByLower 0.001

vtkGeometryFilter geometry
  geometry SetInputConnection [threshold1 GetOutputPort]

vtkImageGradient grad
grad SetDimensionality 3
grad SetInputConnection [sample GetOutputPort]
grad Update

vtkImageMathematics mult
mult SetOperationToMultiply
mult SetInput1Data [sample GetOutput]
mult SetInput2Data [sample GetOutput]

vtkImageToStructuredPoints itosp
itosp SetInputConnection [mult GetOutputPort]
itosp SetVectorInputData [grad GetOutput]
itosp Update

vtkSubPixelPositionEdgels sub
sub SetInputConnection [geometry GetOutputPort]
sub SetGradMapsData [itosp GetOutput]


vtkDataSetMapper mapper
  mapper SetInputConnection [sub GetOutputPort]

vtkActor actor
  actor SetMapper mapper

# Create renderer stuff
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor actor
ren1 ResetCamera
[ren1 GetActiveCamera] Azimuth 20
[ren1 GetActiveCamera] Elevation 30
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 450 450

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}
set cam1 [ren1 GetActiveCamera]
$cam1 Zoom 1.4
iren Initialize


# prevent the tk window from showing up then start the event loop
wm withdraw .
