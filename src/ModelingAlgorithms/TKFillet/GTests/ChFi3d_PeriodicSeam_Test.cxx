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
  // SHA-256: AF03872EB5537ABCE6EA1CC94D867CF15649E139FBAF70C10E4145ADDA80F1A5.
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

TEST(ChFi3d_PeriodicSeamTest, SyntheticSeamSingleAndPairRadiusTwo)
{
  struct FilletCase
  {
    bool useFirst;
    bool useSecond;
    bool reverseOrder;
  };

  const FilletCase aCases[] = {{true, false, false},
                               {false, true, false},
                               {true, true, false},
                               {true, true, true}};
  const double     aRadius  = 20.0;
  for (const FilletCase& aCase : aCases)
  {
    SCOPED_TRACE(::testing::Message() << aCase.useFirst << aCase.useSecond << aCase.reverseOrder);
    TopoDS_Shape aInput;
    BRep_Builder aBuilder;
    ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
    ASSERT_TRUE(BRepCheck_Analyzer(aInput).IsValid());
    const TopoDS_Edge aFirst =
      findEdge(aInput, gp_Pnt(0.0, aRadius, 10.0), gp_Pnt(0.0, aRadius, 20.0));
    const TopoDS_Edge aSecond =
      findEdge(aInput, gp_Pnt(0.0, aRadius, 10.0), gp_Pnt(0.0, 0.0, 10.0));
    ASSERT_FALSE(aFirst.IsNull());
    ASSERT_FALSE(aSecond.IsNull());

    BRepFilletAPI_MakeFillet aFillet(aInput, ChFi3d_Rational);
    if (aCase.reverseOrder)
    {
      aFillet.Add(2.0, aSecond);
      aFillet.Add(2.0, aFirst);
    }
    else
    {
      if (aCase.useFirst)
      {
        aFillet.Add(2.0, aFirst);
      }
      if (aCase.useSecond)
      {
        aFillet.Add(2.0, aSecond);
      }
    }
    aFillet.Build();
    ASSERT_TRUE(aFillet.IsDone());
    ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));
    if (aCase.useFirst)
    {
      checkGeneratedRadius(aFillet, aFirst);
    }
    if (aCase.useSecond)
    {
      checkGeneratedRadius(aFillet, aSecond);
    }
  }
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, RadialPairThenCentralRadiusTwo)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));

  const double      aRadius  = 22.866527557373047;
  const TopoDS_Edge aRadialA = findEdge(aInput, gp_Pnt(0.0, aRadius, 10.0), gp_Pnt(0.0, 0.0, 10.0));
  const TopoDS_Edge aRadialB = findEdge(aInput, gp_Pnt(0.0, 0.0, 10.0), gp_Pnt(aRadius, 0.0, 10.0));
  ASSERT_FALSE(aRadialA.IsNull());
  ASSERT_FALSE(aRadialB.IsNull());

  BRepFilletAPI_MakeFillet aFirst(aInput, ChFi3d_Rational);
  aFirst.Add(2.0, aRadialA);
  aFirst.Add(2.0, aRadialB);
  aFirst.Build();
  ASSERT_TRUE(aFirst.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFirst.Shape()));

  // The first feature consumes the lower 2 mm of the central edge.
  // Resolve its unique successor on the actual first result.
  const TopoDS_Edge aCentral =
    findEdge(aFirst.Shape(), gp_Pnt(0.0, 0.0, 12.0), gp_Pnt(0.0, 0.0, 20.0));
  ASSERT_FALSE(aCentral.IsNull());
  BRepFilletAPI_MakeFillet aSecond(aFirst.Shape(), ChFi3d_Rational);
  aSecond.Add(2.0, aCentral);
  aSecond.Build();
  ASSERT_TRUE(aSecond.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aSecond.Shape()));
  checkGeneratedRadius(aSecond, aCentral);

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(aSecond.Shape(), TopAbs_EDGE, aResultEdges);
  EXPECT_FALSE(aResultEdges.Contains(aCentral));
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, BottomRimDoesNotDuplicateClosedSeam)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  const double      aRadius = 22.866527557373047;
  const TopoDS_Edge aBottom =
    findEdge(aInput, gp_Pnt(0.0, aRadius, 0.0), gp_Pnt(0.0, aRadius, 0.0));
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

TEST(ChFi3d_PeriodicSeamTest, RadialPairAndTripleOrderRadiusTwo)
{
  struct EdgeOrder
  {
    int count;
    int indices[3];
  };

  const EdgeOrder anOrders[] = {{2, {0, 1, 0}},
                                {2, {1, 0, 0}},
                                {3, {0, 1, 2}},
                                {3, {0, 2, 1}},
                                {3, {1, 0, 2}},
                                {3, {1, 2, 0}},
                                {3, {2, 0, 1}},
                                {3, {2, 1, 0}}};
  const double    aRadius    = 22.866527557373047;
  for (const EdgeOrder& anOrder : anOrders)
  {
    SCOPED_TRACE(::testing::Message() << anOrder.count << ": " << anOrder.indices[0] << ","
                                      << anOrder.indices[1] << "," << anOrder.indices[2]);
    TopoDS_Shape aInput;
    BRep_Builder aBuilder;
    ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
    const TopoDS_Edge anEdges[] = {
      findEdge(aInput, gp_Pnt(0.0, aRadius, 10.0), gp_Pnt(0.0, 0.0, 10.0)),
      findEdge(aInput, gp_Pnt(0.0, 0.0, 10.0), gp_Pnt(aRadius, 0.0, 10.0)),
      findEdge(aInput, gp_Pnt(0.0, 0.0, 10.0), gp_Pnt(0.0, 0.0, 20.0))};
    for (const TopoDS_Edge& anEdge : anEdges)
    {
      ASSERT_FALSE(anEdge.IsNull());
    }

    BRepFilletAPI_MakeFillet aFillet(aInput, ChFi3d_Rational);
    for (int anIndex = 0; anIndex < anOrder.count; ++anIndex)
    {
      aFillet.Add(2.0, anEdges[anOrder.indices[anIndex]]);
    }
    aFillet.Build();
    ASSERT_TRUE(aFillet.IsDone());
    ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));
    for (int anIndex = 0; anIndex < anOrder.count; ++anIndex)
    {
      checkGeneratedRadius(aFillet, anEdges[anOrder.indices[anIndex]]);
    }
  }
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, SequentialRadialFeaturesRadiusTwo)
{
  const double aRadius = 22.866527557373047;
  for (int aFirstIndex = 0; aFirstIndex < 2; ++aFirstIndex)
  {
    SCOPED_TRACE(aFirstIndex);
    TopoDS_Shape aInput;
    BRep_Builder aBuilder;
    ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
    const gp_Pnt      aCenter(0.0, 0.0, 10.0);
    const gp_Pnt      anOuter[]  = {gp_Pnt(0.0, aRadius, 10.0), gp_Pnt(aRadius, 0.0, 10.0)};
    const TopoDS_Edge aFirstSeed = findEdge(aInput, aCenter, anOuter[aFirstIndex]);
    ASSERT_FALSE(aFirstSeed.IsNull());

    BRepFilletAPI_MakeFillet aFirst(aInput, ChFi3d_Rational);
    aFirst.Add(2.0, aFirstSeed);
    aFirst.Build();
    ASSERT_TRUE(aFirst.IsDone());
    ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFirst.Shape()));

    // Resolve the opposite radial edge on the first result, not by an
    // index from the original BREP.
    const TopoDS_Edge aSecondSeed = findEdge(aFirst.Shape(), aCenter, anOuter[1 - aFirstIndex]);
    ASSERT_FALSE(aSecondSeed.IsNull());
    BRepFilletAPI_MakeFillet aSecond(aFirst.Shape(), ChFi3d_Rational);
    aSecond.Add(2.0, aSecondSeed);
    aSecond.Build();
    ASSERT_TRUE(aSecond.IsDone());
    ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aSecond.Shape()));
    checkGeneratedRadius(aSecond, aSecondSeed);

    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
    TopExp::MapShapes(aSecond.Shape(), TopAbs_EDGE, aResultEdges);
    EXPECT_FALSE(aResultEdges.Contains(aSecondSeed));
  }
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, UpperStepRimSplitsPeriodicSeamRadiusTwo)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  const double      aRadius = 22.866527557373047;
  const TopoDS_Edge aRim = findEdge(aInput, gp_Pnt(0.0, aRadius, 10.0), gp_Pnt(aRadius, 0.0, 10.0));
  const TopoDS_Edge aSeam = findEdge(aInput, gp_Pnt(0.0, aRadius, 0.0), gp_Pnt(0.0, aRadius, 10.0));
  ASSERT_FALSE(aRim.IsNull());
  ASSERT_FALSE(aSeam.IsNull());

  BRepFilletAPI_MakeFillet aFillet(aInput, ChFi3d_Rational);
  aFillet.Add(2.0, aRim);
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  EXPECT_EQ(aFillet.NbContours(), 1);
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));

  const NCollection_List<TopoDS_Shape>& aModifiedSeam = aFillet.Modified(aSeam);
  ASSERT_FALSE(aModifiedSeam.IsEmpty());
  bool hasLowerSeam = false;
  for (const TopoDS_Shape& aShape : aModifiedSeam)
  {
    if (aShape.ShapeType() == TopAbs_EDGE)
    {
      hasLowerSeam = hasLowerSeam
                     || isEdgeBetween(TopoDS::Edge(aShape),
                                      gp_Pnt(0.0, aRadius, 0.0),
                                      gp_Pnt(0.0, aRadius, 8.0));
    }
  }
  EXPECT_TRUE(hasLowerSeam);

  bool hasRadiusTwoTorus = false;
  for (const TopoDS_Shape& aShape : aFillet.Generated(aRim))
  {
    if (aShape.ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    BRepAdaptor_Surface aSurface(TopoDS::Face(aShape));
    hasRadiusTwoTorus = hasRadiusTwoTorus
                        || (aSurface.GetType() == GeomAbs_Torus
                            && std::abs(aSurface.Torus().MinorRadius() - 2.0) < 1.0e-6);
  }
  EXPECT_TRUE(hasRadiusTwoTorus);
}

//==================================================================================================

TEST(ChFi3d_PeriodicSeamTest, AllPhysicalEdgesAcrossPeriodicCornersRadiusTwo)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  const double aRadius         = 22.866527557373047;
  const gp_Pnt aEndpoints[][2] = {{{0.0, aRadius, 0.0}, {0.0, aRadius, 0.0}},
                                  {{0.0, aRadius, 10.0}, {aRadius, 0.0, 10.0}},
                                  {{aRadius, 0.0, 10.0}, {aRadius, 0.0, 20.0}},
                                  {{aRadius, 0.0, 20.0}, {0.0, aRadius, 20.0}},
                                  {{0.0, aRadius, 10.0}, {0.0, aRadius, 20.0}},
                                  {{0.0, aRadius, 10.0}, {0.0, 0.0, 10.0}},
                                  {{0.0, 0.0, 10.0}, {aRadius, 0.0, 10.0}},
                                  {{0.0, 0.0, 10.0}, {0.0, 0.0, 20.0}},
                                  {{0.0, 0.0, 20.0}, {aRadius, 0.0, 20.0}},
                                  {{0.0, aRadius, 20.0}, {0.0, 0.0, 20.0}}};

  BRepFilletAPI_MakeFillet aFillet(aInput, ChFi3d_Rational);
  TopoDS_Edge              aSeeds[10];
  for (int anIndex = 0; anIndex < 10; ++anIndex)
  {
    aSeeds[anIndex] = findEdge(aInput, aEndpoints[anIndex][0], aEndpoints[anIndex][1]);
    ASSERT_FALSE(aSeeds[anIndex].IsNull()) << anIndex;
    aFillet.Add(2.0, aSeeds[anIndex]);
  }
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  ASSERT_NO_FATAL_FAILURE(checkCompleteSolid(aFillet.Shape()));

  for (const TopoDS_Edge& aSeed : aSeeds)
  {
    const NCollection_List<TopoDS_Shape>& aGenerated = aFillet.Generated(aSeed);
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
                         && std::abs(aSurface.Cylinder().Radius() - 2.0) < 1.0e-6)
                     || (aSurface.GetType() == GeomAbs_Torus
                         && std::abs(aSurface.Torus().MinorRadius() - 2.0) < 1.0e-6);
    }
    EXPECT_TRUE(hasRadiusTwo);
  }

  int aSplineCorners = 0;
  for (TopExp_Explorer aFaceIt(aFillet.Shape(), TopAbs_FACE); aFaceIt.More(); aFaceIt.Next())
  {
    const TopoDS_Face   aFace = TopoDS::Face(aFaceIt.Current());
    BRepAdaptor_Surface aSurface(aFace);
    if (aSurface.GetType() != GeomAbs_BSplineSurface)
    {
      continue;
    }
    ++aSplineCorners;
    double uMin, uMax, vMin, vMax;
    BRepTools::UVBounds(aFace, uMin, uMax, vMin, vMax);
    gp_Pnt aPoint;
    gp_Vec aDU, aDV;
    aSurface.D1((uMin + uMax) / 2.0, (vMin + vMax) / 2.0, aPoint, aDU, aDV);
    const gp_Vec aReferenceNormal = aDU.Crossed(aDV);
    ASSERT_GT(aReferenceNormal.Magnitude(), 1.0e-2);
    for (int aUIndex = 1; aUIndex < 20; ++aUIndex)
    {
      for (int aVIndex = 1; aVIndex < 20; ++aVIndex)
      {
        const double u = uMin + (uMax - uMin) * aUIndex / 20.0;
        const double v = vMin + (vMax - vMin) * aVIndex / 20.0;
        aSurface.D1(u, v, aPoint, aDU, aDV);
        const gp_Vec aNormal = aDU.Crossed(aDV);
        EXPECT_GT(aNormal.Magnitude(), 1.0e-2);
        EXPECT_GT(aNormal.Dot(aReferenceNormal), 0.0);
      }
    }
  }
  EXPECT_EQ(aSplineCorners, 2);
}
