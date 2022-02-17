package require vtk
package require vtkinteraction

# Demonstrate how to use the vtkInteractorEventRecorder to play back some
# events.

# Create a mace out of filters.
#
vtkSphereSource sphere
vtkConeSource cone
vtkGlyph3D glyph
    glyph SetInputConnection [sphere GetOutputPort]
    glyph SetSourceConnection [cone GetOutputPort]
    glyph SetVectorModeToUseNormal
    glyph SetScaleModeToScaleByVector
    glyph SetScaleFactor 0.25

# The sphere and spikes are appended into a single polydata. This just makes things
# simpler to manage.
vtkAppendPolyData apd
    apd AddInputConnection [glyph GetOutputPort]
    apd AddInputConnection [sphere GetOutputPort]
vtkPolyDataMapper maceMapper
    maceMapper SetInputConnection [apd GetOutputPort]
vtkLODActor maceActor
    maceActor SetMapper maceMapper
    maceActor VisibilityOn

# This portion of the code clips the mace with the vtkPlanes implicit function.
# The clipped region is colored green.
vtkPlanes planes
vtkClipPolyData clipper
    clipper SetInputConnection [apd GetOutputPort]
    clipper SetClipFunction planes
    clipper InsideOutOn
vtkPolyDataMapper selectMapper
    selectMapper SetInputConnection [clipper GetOutputPort]
vtkLODActor selectActor
    selectActor SetMapper selectMapper
    [selectActor GetProperty] SetColor 0 1 0
    selectActor VisibilityOff
    selectActor SetScale 1.01 1.01 1.01

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1

vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin
    iren AddObserver ExitEvent {exit}

# The SetInteractor method is how 3D widgets are associated with the render
# window interactor. Internally, SetInteractor sets up a bunch of callbacks
# using the Command/Observer mechanism (AddObserver()).
ren1 AddActor maceActor
ren1 AddActor selectActor

# Add the actors to the renderer, set the background and size
#
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 300 300

# Place the interactor initially. The input to a 3D widget is used to
# initially position and scale the widget. The EndInteractionEvent is
# observed which invokes the SelectPolygons callback.
vtkBoxRepresentation boxRep
boxRep SetPlaceFactor  0.75
eval boxRep PlaceWidget [[glyph GetOutput] GetBounds]

vtkBoxWidget2  boxWidget
boxWidget SetInteractor  iren
boxWidget SetRepresentation  boxRep
boxWidget AddObserver EndInteractionEvent SelectPolygons
boxWidget SetPriority 1

# record events
vtkInteractorEventRecorder recorder
recorder SetInteractor iren
recorder SetFileName "$VTK_DATA_ROOT/Data/EventRecording.log"

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize
renWin Render
#recorder Record

# prevent the tk window from showing up then start the event loop
wm withdraw .

# This does the actual work: updates the vtkPlanes implicit function.
# This in turn causes the pipeline to update.
proc SelectPolygons {} {
   boxRep GetPlanes planes
   selectActor VisibilityOn
}

recorder Play
recorder Off
set threshold 15

