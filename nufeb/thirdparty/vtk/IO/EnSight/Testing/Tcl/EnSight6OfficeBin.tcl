package require vtk
package require vtkinteraction

vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# read data
#
vtkGenericEnSightReader reader
# Make sure all algorithms use the composite data pipeline
vtkCompositeDataPipeline cdp
reader SetDefaultExecutivePrototype cdp
cdp Delete
reader SetCaseFileName "$VTK_DATA_ROOT/Data/EnSight/office6_bin.case"
reader Update

vtkStructuredGridOutlineFilter outline
    outline SetInputConnection [reader GetOutputPort]
vtkHierarchicalPolyDataMapper mapOutline
    mapOutline SetInputConnection [outline GetOutputPort]
vtkActor outlineActor
    outlineActor SetMapper mapOutline
    [outlineActor GetProperty] SetColor 0 0 0

# Create source for streamtubes
vtkStreamTracer streamer
    streamer SetInputConnection [reader GetOutputPort]
    streamer SetStartPosition 0.1 2.1 0.5
    streamer SetMaximumPropagation 500
    streamer SetInitialIntegrationStep 0.1
    streamer SetIntegrationDirectionToForward

vtkConeSource cone
    cone SetResolution 8
vtkGlyph3D cones
    cones SetInputConnection [streamer GetOutputPort]
    cones SetSourceConnection [cone GetOutputPort]
    cones SetScaleFactor 3
    cones SetInputArrayToProcess 1 0 0 0 "vectors"
    cones SetScaleModeToScaleByVector

vtkHierarchicalPolyDataMapper mapCones
    mapCones SetInputConnection [cones GetOutputPort]
eval mapCones SetScalarRange [[[reader GetOutput] GetBlock 0] GetScalarRange]
vtkActor conesActor
    conesActor SetMapper mapCones

ren1 AddActor outlineActor
ren1 AddActor conesActor

ren1 SetBackground 0.4 0.4 0.5

renWin SetSize 300 300
iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize

# interact with data
wm withdraw .

reader SetDefaultExecutivePrototype {}

