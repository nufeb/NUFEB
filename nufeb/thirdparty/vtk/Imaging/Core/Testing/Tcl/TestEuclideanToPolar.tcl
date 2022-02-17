# This Script test the euclidean to polar by coverting 2D vectors
# from a gradient into polar, which is converted into HSV, and then to RGB.
package require vtk




# Image pipeline

vtkImageGaussianSource gauss
gauss SetWholeExtent 0 255 0 255 0 44
gauss SetCenter 127 127 22
gauss SetStandardDeviation 50.0
gauss SetMaximum 10000.0

vtkImageGradient gradient
gradient SetInputConnection [gauss GetOutputPort]
gradient SetDimensionality 2

vtkImageEuclideanToPolar polar
polar SetInputConnection [gradient GetOutputPort]
polar SetThetaMaximum 255

vtkImageViewer viewer
#viewer DebugOn
viewer SetInputConnection [polar GetOutputPort]
viewer SetZSlice 22
viewer SetColorWindow 255
viewer SetColorLevel 127.5


#make interface
viewer Render







