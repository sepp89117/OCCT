// Copyright (c) 2026 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

// Regression for OCC28354: fillet of a horizontal edge touching a cylindrical seam.
// Fixture: exact frozen official bug28354.brep, SHA-256
// 847AB883D542419D01C68E3A45B482B5E3A75D6B2ECBBDD05E957AC15C781253.
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepAlgoAPI_Check.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_Shape.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

namespace
{
TCollection_AsciiString fixturePath()
{
  TCollection_AsciiString aPath(__FILE__);
  const int aSlash = std::max(aPath.SearchFromEnd("/"), aPath.SearchFromEnd("\\"));
  aPath.Trunc(aSlash);
  aPath += "data/bug28354.brep";
  return aPath;
}

TopoDS_Edge findHorizontalEdgeAtSeam(const TopoDS_Shape& shape)
{
  const gp_Pnt a(-3.46410161513775, 2.0, 4.0);
  const gp_Pnt b( 3.46410161513775, 2.0, 4.0);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
  TopExp::MapShapes(shape, TopAbs_EDGE, edges);
  TopoDS_Edge found;
  for (int i = 1; i <= edges.Extent(); ++i)
  {
    const TopoDS_Edge edge = TopoDS::Edge(edges(i));
    TopoDS_Vertex first, last;
    TopExp::Vertices(edge, first, last);
    if (first.IsNull() || last.IsNull())
      continue;
    const gp_Pnt p = BRep_Tool::Pnt(first);
    const gp_Pnt q = BRep_Tool::Pnt(last);
    const bool match = (p.Distance(a) < 1.0e-6 && q.Distance(b) < 1.0e-6)
                       || (p.Distance(b) < 1.0e-6 && q.Distance(a) < 1.0e-6);
    if (match)
    {
      if (!found.IsNull())
        return TopoDS_Edge();
      found = edge;
    }
  }
  return found;
}
}

TEST(ChFi3dExternalBadPaths, OCC28354SeamFilletProducesValidSolid)
{
  TopoDS_Shape input;
  BRep_Builder builder;
  ASSERT_TRUE(BRepTools::Read(input, fixturePath().ToCString(), builder));
  ASSERT_TRUE(BRepCheck_Analyzer(input).IsValid());
  const TopoDS_Edge seed = findHorizontalEdgeAtSeam(input);
  ASSERT_FALSE(seed.IsNull());

  BRepFilletAPI_MakeFillet fillet(input);
  fillet.SetParams(0.01, 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-5, 1.0e-3);
  fillet.SetContinuity(GeomAbs_C1, 0.01);
  fillet.Add(2.2, seed); // The official script's second, formerly BAD construction.
  EXPECT_NO_THROW(fillet.Build());
  ASSERT_TRUE(fillet.IsDone());

  const TopoDS_Shape& result = fillet.Shape();
  ASSERT_FALSE(result.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(result).IsValid());
  BRepAlgoAPI_Check booleanCheck(result, true, true);
  ASSERT_TRUE(booleanCheck.IsValid());
  int solidCount = 0;
  for (TopExp_Explorer it(result, TopAbs_SOLID); it.More(); it.Next())
    ++solidCount;
  EXPECT_EQ(solidCount, 1);

  GProp_GProps inputMass, resultMass;
  BRepGProp::VolumeProperties(input, inputMass);
  BRepGProp::VolumeProperties(result, resultMass);
  EXPECT_GT(resultMass.Mass(), 0.0);
  EXPECT_GT(std::abs(resultMass.Mass() - inputMass.Mass()), 1.0e-5);

  Bnd_Box bounds;
  BRepBndLib::Add(result, bounds);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  EXPECT_NEAR(zmin, 0.0, 0.01); // The fillet at z=4 must retain the cylinder ends.
  EXPECT_NEAR(zmax, 8.0, 0.01);
}
