package require vtk
package require vtkinteraction

vtkPoints tetraPoints
  tetraPoints SetNumberOfPoints 4
  tetraPoints InsertPoint 0 0 0 0
  tetraPoints InsertPoint 1 1 0 0
  tetraPoints InsertPoint 2 .5 1 0
  tetraPoints InsertPoint 3 .5 .5 1

vtkTetra aTetra
  [aTetra GetPointIds] SetId 0 0
  [aTetra GetPointIds] SetId 1 1
  [aTetra GetPointIds] SetId 2 2
  [aTetra GetPointIds] SetId 3 3

vtkUnstructuredGrid aTetraGrid
  aTetraGrid Allocate 1 1
  aTetraGrid InsertNextCell [aTetra GetCellType] [aTetra GetPointIds]
  aTetraGrid SetPoints tetraPoints

vtkSubdivideTetra sub
  sub SetInputData aTetraGrid

vtkShrinkFilter shrinker
  shrinker SetInputConnection [sub GetOutputPort]

vtkDataSetMapper mapper
  mapper SetInputConnection [shrinker GetOutputPort]

vtkActor actor
  actor SetMapper mapper
  [actor GetProperty] SetColor 0.7400 0.9900 0.7900

# define graphics stuff
vtkRenderer ren1
vtkRenderWindow renWin
  renWin AddRenderer ren1
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin

ren1 AddActor actor
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 300 300

set cam1 [ren1 GetActiveCamera]
$cam1 SetClippingRange 0.183196 9.15979
$cam1 SetFocalPoint 0.579471 0.462507 0.283392
$cam1 SetPosition -1.04453 0.345281 -0.556222
$cam1 SetViewUp 0.197321 0.843578 -0.499441
ren1 ResetCameraClippingRange

renWin Render

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize
wm withdraw .



