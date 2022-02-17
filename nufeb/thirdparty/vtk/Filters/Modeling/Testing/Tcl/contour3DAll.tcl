package require vtk
package require vtkinteraction
package require vtktesting

# On older Macs, 10 is too low. Due to what looks like a driver bug
# spectral lighting behaves sort of weird and produces small differences
set threshold 30

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create pipeline
#
vtkStructuredPointsReader slc
  slc SetFileName "$VTK_DATA_ROOT/Data/ironProt.vtk"

set colors "$flesh $banana $grey $pink $carrot $gainsboro $tomato $gold $thistle $chocolate"
set types  "UnsignedChar Char Short UnsignedShort Int UnsignedInt Long UnsignedLong Float Double"
set i 1
set c 0
foreach vtkType $types {
  vtkImageClip clip$vtkType
    clip$vtkType SetInputConnection [slc GetOutputPort]
    clip$vtkType SetOutputWholeExtent -1000 1000 -1000 1000 $i [expr $i + 5]
  incr i 5
  vtkImageCast castTo$vtkType
    castTo$vtkType SetOutputScalarTypeTo$vtkType
    castTo$vtkType SetInputConnection [clip$vtkType GetOutputPort]
    castTo$vtkType ClampOverflowOn

  vtkContourFilter iso$vtkType
    iso$vtkType SetInputConnection [castTo$vtkType GetOutputPort]
    iso$vtkType GenerateValues 1 30 30
    iso$vtkType ComputeScalarsOff
    iso$vtkType ComputeGradientsOff

  vtkPolyDataMapper iso${vtkType}Mapper
    iso${vtkType}Mapper SetInputConnection [iso$vtkType GetOutputPort]
    iso${vtkType}Mapper ImmediateModeRenderingOn

  vtkActor iso${vtkType}Actor
    iso${vtkType}Actor SetMapper iso${vtkType}Mapper
    [iso${vtkType}Actor GetProperty] SetDiffuseColor [lindex $colors $c] [lindex $colors [expr $c + 1]] [lindex $colors [expr $c + 2]]
    [iso${vtkType}Actor GetProperty] SetSpecularPower 30
    [iso${vtkType}Actor GetProperty] SetDiffuse .7
    [iso${vtkType}Actor GetProperty] SetSpecular .5
    incr c 3
    ren1 AddActor iso${vtkType}Actor
}

vtkOutlineFilter outline
  outline SetInputConnection [slc GetOutputPort]
vtkPolyDataMapper outlineMapper
  outlineMapper SetInputConnection [outline GetOutputPort]
vtkActor outlineActor
  outlineActor SetMapper outlineMapper
  outlineActor VisibilityOff

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 SetBackground 0.9 .9 .9
ren1 ResetCamera
[ren1 GetActiveCamera] SetViewAngle 30
[ren1 GetActiveCamera] Elevation 20
[ren1 GetActiveCamera] Azimuth 20
[ren1 GetActiveCamera] Zoom 1.5
ren1 ResetCameraClippingRange

renWin SetSize 400 400
iren Initialize

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .
