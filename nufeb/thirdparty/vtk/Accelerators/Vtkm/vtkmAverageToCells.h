//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2012 Sandia Corporation.
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//=============================================================================

#ifndef vtkmAverageToCells_h
#define vtkmAverageToCells_h

#include "vtkAcceleratorsVTKmModule.h" //required for correct implementation
#include "vtkDataSetAlgorithm.h"

class VTKACCELERATORSVTKM_EXPORT vtkmAverageToCells : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkmAverageToCells, vtkDataSetAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;
  static vtkmAverageToCells* New();

protected:
  vtkmAverageToCells();
  ~vtkmAverageToCells();

  virtual int RequestData(vtkInformation*, vtkInformationVector**,
                          vtkInformationVector*) VTK_OVERRIDE;

private:
  vtkmAverageToCells(const vtkmAverageToCells&) VTK_DELETE_FUNCTION;
  void operator=(const vtkmAverageToCells&) VTK_DELETE_FUNCTION;
};

#endif // vtkmAverageToCells_h
// VTK-HeaderTest-Exclude: vtkmAverageToCells.h
