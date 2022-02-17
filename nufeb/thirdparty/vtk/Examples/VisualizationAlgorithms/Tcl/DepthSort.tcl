# This example demonstrates the use of vtkDepthSortPolyData. This is a
# poor man's algorithm to sort polygons for proper transparent blending.
# It sorts polygons based on a single point (i.e., centroid) so the sorting
# may not work for overlapping or intersection polygons.
#
package require vtk
package require vtkinteraction

# Create a bunch of spheres that overlap and cannot be easily arranged
# so that the blending works without sorting. They are appended into a
# single vtkPolyData because the filter only sorts within a single
# vtkPolyData input.
#
vtkSphereSource sphere
    sphere SetThetaResolution 80
    sphere SetPhiResolution 40
    sphere SetRadius 1
    sphere SetCenter 0 0 0
vtkSphereSource sphere2
    sphere2 SetThetaResolution 80
    sphere2 SetPhiResolution 40
    sphere2 SetRadius 0.5
    sphere2 SetCenter 1 0 0
vtkSphereSource sphere3
    sphere3 SetThetaResolution 80
    sphere3 SetPhiResolution 40
    sphere3 SetRadius 0.5
    sphere3 SetCenter -1 0 0
vtkSphereSource sphere4
    sphere4 SetThetaResolution 80
    sphere4 SetPhiResolution 40
    sphere4 SetRadius 0.5
    sphere4 SetCenter 0 1 0
vtkSphereSource sphere5
    sphere5 SetThetaResolution 80
    sphere5 SetPhiResolution 40
    sphere5 SetRadius 0.5
    sphere5 SetCenter 0 -1 0
vtkAppendPolyData appendData
    appendData AddInputConnection [sphere GetOutputPort]
    appendData AddInputConnection [sphere2 GetOutputPort]
    appendData AddInputConnection [sphere3 GetOutputPort]
    appendData AddInputConnection [sphere4 GetOutputPort]
    appendData AddInputConnection [sphere5 GetOutputPort]

# The dephSort object is set up to generate scalars representing
# the sort depth.  A camera is assigned for the sorting. The camera
# define the sort vector (position and focal point).
vtkCamera camera
vtkDepthSortPolyData depthSort
    depthSort SetInputConnection [appendData GetOutputPort]
    depthSort SetDirectionToBackToFront
    depthSort SetVector 1 1 1
    depthSort SetCamera camera
    depthSort SortScalarsOn
    depthSort Update

vtkPolyDataMapper mapper
    mapper SetInputConnection [depthSort GetOutputPort]
    mapper SetScalarRange 0 [[depthSort GetOutput] GetNumberOfCells]
vtkActor actor
    actor SetMapper mapper
    [actor GetProperty] SetOpacity 0.5
    [actor GetProperty] SetColor 1 0 0
    actor RotateX -72

# If an Prop3D is supplied, then its transformation matrix is taken
# into account during sorting.
depthSort SetProp3D actor

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
    ren1 SetActiveCamera camera
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor actor
ren1 SetBackground 1 1 1
renWin SetSize 300 200

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}
ren1 ResetCamera
[ren1 GetActiveCamera] Zoom 2.2
renWin Render

# prevent the tk window from showing up then start the event loop
wm withdraw .

iren Start