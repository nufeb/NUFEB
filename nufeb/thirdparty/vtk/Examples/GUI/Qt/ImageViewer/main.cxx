

/*=========================================================================
  Copyright 2004 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/

/*========================================================================
 For general information about using VTK and Qt, see:
 http://www.trolltech.com/products/3rdparty/vtksupport.html
=========================================================================*/

/*========================================================================
 !!! WARNING for those who want to contribute code to this file.
 !!! If you use a commercial edition of Qt, you can modify this code.
 !!! If you use an open source version of Qt, you are free to modify
 !!! and use this code within the guidelines of the GPL license.
 !!! Unfortunately, you cannot contribute the changes back into this
 !!! file.  Doing so creates a conflict between the GPL and BSD-like VTK
 !!! license.
=========================================================================*/

#include <QApplication>

#include "vtkImageViewer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkPNGReader.h"
#include "vtkTestUtilities.h"

#include "QVTKWidget.h"

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

  QVTKWidget widget;
  widget.resize(256,256);

  vtkPNGReader* reader = vtkPNGReader::New();
  char* fname = vtkTestUtilities::ExpandDataFileName(argc, argv, "Data/vtk.png");
  reader->SetFileName(fname);
  delete [] fname;

  vtkImageViewer* image_view = vtkImageViewer::New();
  image_view->SetInputConnection(reader->GetOutputPort());

  widget.SetRenderWindow(image_view->GetRenderWindow());
  image_view->SetupInteractor(widget.GetRenderWindow()->GetInteractor());

  image_view->SetColorLevel(138.5);
  image_view->SetColorWindow(233);

  widget.show();

  app.exec();


  image_view->Delete();
  reader->Delete();

  return 0;
}


