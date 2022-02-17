package require vtk

# Shift and scale an image (in that order)
# This filter is useful for converting to a lower precision data type.

vtkImageReader reader
[reader GetExecutive] SetReleaseDataFlag 0 0
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 63 0 63 1 93
reader SetFilePrefix "$VTK_DATA_ROOT/Data/headsq/quarter"
reader SetDataMask 0x7fff

vtkImageShiftScale shiftScale
shiftScale SetInputConnection [reader GetOutputPort]
shiftScale SetShift 0
shiftScale SetScale 0.5
shiftScale SetOutputScalarTypeToDouble

vtkImageShiftScale shiftScale2
shiftScale2 SetInputConnection [shiftScale GetOutputPort]
shiftScale2 SetShift 0
shiftScale2 SetScale 2.0

vtkImageMagnify mag
  mag SetInputConnection [shiftScale2 GetOutputPort]
  mag SetMagnificationFactors 4 4 1
  mag InterpolateOff

vtkImageViewer viewer
viewer SetInputConnection [mag GetOutputPort]
viewer SetColorWindow 1024
viewer SetColorLevel 1024

#make interface
source [file join [file dirname [info script]] WindowLevelInterface.tcl]


#vtkPNMWriter w
#w SetFileName "D:/vtknew/vtk/graphics/examplesTcl/mace2.ppm"
#w SetInputConnection [shiftScale GetOutputPort]
#w Write
