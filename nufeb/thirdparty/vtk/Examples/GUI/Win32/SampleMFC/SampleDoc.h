/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SampleDoc.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// SampleDoc.h : interface of the CSampleDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_SAMPLEDOC_H__B7F7B861_EEC9_11D2_87FE_0060082B79FD__INCLUDED_)
#define AFX_SAMPLEDOC_H__B7F7B861_EEC9_11D2_87FE_0060082B79FD__INCLUDED_

#ifdef _MSC_VER
#pragma once
#endif

#include "vtkMFCDocument.h"
#include "vtkDataSetReader.h"
#include "vtkDataSetMapper.h"

class CSampleDoc : public vtkMFCDocument
{
protected: // create from serialization only
        CSampleDoc();
        DECLARE_DYNCREATE(CSampleDoc)

// Attributes
public:

// Operations
public:

// Overrides
        // ClassWizard generated virtual function overrides
        //{{AFX_VIRTUAL(CSampleDoc)
        public:
        virtual BOOL OnNewDocument();
        virtual void Serialize(CArchive& ar);
        virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
        //}}AFX_VIRTUAL

// Implementation
public:
        virtual ~CSampleDoc();
#ifdef _DEBUG
        virtual void AssertValid() const;
        virtual void Dump(CDumpContext& dc) const;
#endif

protected:
  vtkDataSetReader *Reader;
  vtkDataSetMapper *Mapper;
  vtkActor *Actor;

// Generated message map functions
protected:
        //{{AFX_MSG(CSampleDoc)
                // NOTE - the ClassWizard will add and remove member functions here.
                //    DO NOT EDIT what you see in these blocks of generated code !
        //}}AFX_MSG
        DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SAMPLEDOC_H__B7F7B861_EEC9_11D2_87FE_0060082B79FD__INCLUDED_)
