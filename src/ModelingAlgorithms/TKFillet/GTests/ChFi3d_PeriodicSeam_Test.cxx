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

#include <Bnd_Box.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Check.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <ChFi3d_FilletShape.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <gp_Vec.hxx>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

namespace
{
//==================================================================================================

bool isEdgeBetween(const TopoDS_Edge& theEdge, const gp_Pnt& theFirst, const gp_Pnt& theLast)
{
  TopoDS_Vertex aFirstVertex;
  TopoDS_Vertex aLastVertex;
  TopExp::Vertices(theEdge, aFirstVertex, aLastVertex);
  if (aFirstVertex.IsNull() || aLastVertex.IsNull())
  {
    return false;
  }
  const gp_Pnt aFirst = BRep_Tool::Pnt(aFirstVertex);
  const gp_Pnt aLast  = BRep_Tool::Pnt(aLastVertex);
  const double aTol   = 1.0e-6;
  return (aFirst.Distance(theFirst) < aTol && aLast.Distance(theLast) < aTol)
         || (aFirst.Distance(theLast) < aTol && aLast.Distance(theFirst) < aTol);
}

//==================================================================================================

TopoDS_Edge findEdge(const TopoDS_Shape& theShape, const gp_Pnt& theFirst, const gp_Pnt& theLast)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anEdges;
  TopExp::MapShapes(theShape, TopAbs_EDGE, anEdges);
  TopoDS_Edge aFound;
  for (int anIndex = 1; anIndex <= anEdges.Extent(); ++anIndex)
  {
    const TopoDS_Edge anEdge = TopoDS::Edge(anEdges(anIndex));
    if (isEdgeBetween(anEdge, theFirst, theLast))
    {
      if (!aFound.IsNull())
      {
        return TopoDS_Edge();
      }
      aFound = anEdge;
    }
  }
  return aFound;
}

//==================================================================================================

TCollection_AsciiString fixturePath()
{
  // Public synthetic cylinder and quarter-boss with a real periodic seam.
  TCollection_AsciiString aPath(__FILE__);
  const int               aSlash = std::max(aPath.SearchFromEnd("/"), aPath.SearchFromEnd("\\"));
  aPath.Trunc(aSlash);
  aPath += "data/PeriodicSeamSynthetic.brep";
  return aPath;
}

//==================================================================================================

void checkGeneratedRadius(BRepFilletAPI_MakeFillet& theFillet, const TopoDS_Edge& theSeed)
{
  const NCollection_List<TopoDS_Shape>& aGenerated = theFillet.Generated(theSeed);
  ASSERT_FALSE(aGenerated.IsEmpty());
  bool hasRadiusTwo = false;
  for (const TopoDS_Shape& aShape : aGenerated)
  {
    if (aShape.ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    BRepAdaptor_Surface aSurface(TopoDS::Face(aShape));
    hasRadiusTwo = hasRadiusTwo
                   || (aSurface.GetType() == GeomAbs_Cylinder
                       && std::abs(aSurface.Cylinder().Radius() - 2.0) < 1.0e-6);
  }
  EXPECT_TRUE(hasRadiusTwo);
}

//==================================================================================================

void checkCompleteSolid(const TopoDS_Shape& theShape)
{
  ASSERT_FALSE(theShape.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(theShape).IsValid());
  BRepAlgoAPI_Check aCheck(theShape, true, true);
  ASSERT_TRUE(aCheck.IsValid());

  int aSolids = 0;
  for (TopExp_Explorer aSolidIt(theShape, TopAbs_SOLID); aSolidIt.More(); aSolidIt.Next())
  {
    ++aSolids;
  }
  ASSERT_EQ(aSolids, 1);

  Bnd_Box aBounds;
  BRepBndLib::AddOptimal(theShape, aBounds, false, false);
  double aXMin, aYMin, aZMin, aXMax, aYMax, aZMax;
  aBounds.Get(aXMin, aYMin, aZMin, aXMax, aYMax, aZMax);
  EXPECT_NEAR(aZMin, 0.0, 1.0e-5);
  EXPECT_NEAR(aZMax, 20.0, 1.0e-5);

  GProp_GProps aProperties;
  BRepGProp::VolumeProperties(theShape, aProperties);
  EXPECT_GT(aProperties.Mass(), 15000.0);
  EXPECT_LT(aProperties.Mass(), 16000.0);

  BRepClass3d_SolidClassifier aClassifier(theShape);
  aClassifier.Perform(gp_Pnt(0.0, 0.0, 5.0), 1.0e-7);
  EXPECT_EQ(aClassifier.State(), TopAbs_IN);
  aClassifier.Perform(gp_Pnt(10.0, 10.0, 5.0), 1.0e-7);
  EXPECT_EQ(aClassifier.State(), TopAbs_IN);

  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEdgeFaces;
  TopExp::MapShapesAndAncestors(theShape, TopAbs_EDGE, TopAbs_FACE, anEdgeFaces);
  for (int anIndex = 1; anIndex <= anEdgeFaces.Extent(); ++anIndex)
  {
    const TopoDS_Edge anEdge = TopoDS::Edge(anEdgeFaces.FindKey(anIndex));
    if (!BRep_Tool::Degenerated(anEdge))
    {
      EXPECT_EQ(anEdgeFaces(anIndex).Size(), 2u);
    }
    EXPECT_LT(BRep_Tool::Tolerance(anEdge), 1.0e-3);
  }
}
} // namespace

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, BottomRimDoesNotDuplicateClosedSeam)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  ASSERT_TRUE(BRepCheck_Analyzer(aInput).IsValid());
  const TopoDS_Edge aBottom =
    findEdge(aInput, gp_Pnt(0.0, 20.0, 0.0), gp_Pnt(0.0, 20.0, 0.0));
  ASSERT_FALSE(aBottom.IsNull());

  BRepFilletAPI_MakeFillet aFillet(aInput, ChFi3d_Rational);
  aFillet.Add(2.0, aBottom);
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(aFillet.Shape(), TopAbs_EDGE, aResultEdges);
  EXPECT_FALSE(aResultEdges.Contains(aBottom));
  bool hasRadiusTwo = false;
  for (const TopoDS_Shape& aShape : aFillet.Generated(aBottom))
  {
    if (aShape.ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    BRepAdaptor_Surface aSurface(TopoDS::Face(aShape));
    hasRadiusTwo = hasRadiusTwo
                   || (aSurface.GetType() == GeomAbs_Torus
                       && std::abs(aSurface.Torus().MinorRadius() - 2.0) < 1.0e-6);
  }
  EXPECT_TRUE(hasRadiusTwo);
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, UpperQuarterRimAtPeriodicSeamRadiusTwo)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  ASSERT_TRUE(BRepCheck_Analyzer(aInput).IsValid());
  const TopoDS_Edge anUpperRim =
    findEdge(aInput, gp_Pnt(20.0, 0.0, 20.0), gp_Pnt(0.0, 20.0, 20.0));
  ASSERT_FALSE(anUpperRim.IsNull());

  BRepFilletAPI_MakeFillet aFillet(aInput);
  aFillet.Add(2.0, anUpperRim);
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));

  bool hasRadiusTwo = false;
  for (const TopoDS_Shape& aShape : aFillet.Generated(anUpperRim))
  {
    if (aShape.ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    BRepAdaptor_Surface aSurface(TopoDS::Face(aShape));
    hasRadiusTwo = hasRadiusTwo
                   || (aSurface.GetType() == GeomAbs_Torus
                       && std::abs(aSurface.Torus().MinorRadius() - 2.0) < 1.0e-6);
  }
  EXPECT_TRUE(hasRadiusTwo);

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(aFillet.Shape(), TopAbs_EDGE, aResultEdges);
  EXPECT_FALSE(aResultEdges.Contains(anUpperRim));
}

//==================================================================================================
TEST(ChFi3d_PeriodicSeamTest, UpperRimAndRadialAtPeriodicSeamRadiusTwo)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  ASSERT_TRUE(BRepCheck_Analyzer(aInput).IsValid());
  const TopoDS_Edge anUpperRim =
    findEdge(aInput, gp_Pnt(20.0, 0.0, 20.0), gp_Pnt(0.0, 20.0, 20.0));
  const TopoDS_Edge aRadial =
    findEdge(aInput, gp_Pnt(0.0, 0.0, 20.0), gp_Pnt(0.0, 20.0, 20.0));
  ASSERT_FALSE(anUpperRim.IsNull());
  ASSERT_FALSE(aRadial.IsNull());

  BRepFilletAPI_MakeFillet aFillet(aInput);
  aFillet.Add(2.0, anUpperRim);
  aFillet.Add(2.0, aRadial);
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));
  checkGeneratedRadius(aFillet, aRadial);
  ASSERT_FALSE(aFillet.Generated(anUpperRim).IsEmpty());

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(aFillet.Shape(), TopAbs_EDGE, aResultEdges);
  EXPECT_FALSE(aResultEdges.Contains(anUpperRim));
  EXPECT_FALSE(aResultEdges.Contains(aRadial));
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, LowerAndUpperRimsPreserveCompleteSolidRadiusTwo)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  ASSERT_TRUE(BRepCheck_Analyzer(aInput).IsValid());
  const TopoDS_Edge aLowerRim =
    findEdge(aInput, gp_Pnt(0.0, 20.0, 0.0), gp_Pnt(0.0, 20.0, 0.0));
  const TopoDS_Edge anUpperRim =
    findEdge(aInput, gp_Pnt(20.0, 0.0, 20.0), gp_Pnt(0.0, 20.0, 20.0));
  ASSERT_FALSE(aLowerRim.IsNull());
  ASSERT_FALSE(anUpperRim.IsNull());

  BRepFilletAPI_MakeFillet aFillet(aInput);
  aFillet.Add(2.0, aLowerRim);
  aFillet.Add(2.0, anUpperRim);
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));
  ASSERT_FALSE(aFillet.Generated(aLowerRim).IsEmpty());
  ASSERT_FALSE(aFillet.Generated(anUpperRim).IsEmpty());

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(aFillet.Shape(), TopAbs_EDGE, aResultEdges);
  EXPECT_FALSE(aResultEdges.Contains(aLowerRim));
  EXPECT_FALSE(aResultEdges.Contains(anUpperRim));
}
