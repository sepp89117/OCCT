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
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace
{
// The unmodified 431,177-byte BREP is from OCCT issue #1761, attachment
// conge-y-nervure-base.brep.txt, retrieved 2026-11-21.
// SHA-256: 985CB51C2395CD7C67E670796A8FCE3411B092E5CF75E2492B1B7F983FB691C2.
TCollection_AsciiString fixturePath()
{
  TCollection_AsciiString aPath(__FILE__);
  const int               aSlash = std::max(aPath.SearchFromEnd("/"), aPath.SearchFromEnd("\\"));
  aPath.Trunc(aSlash);
  aPath += "data/Issue1761_YRib.brep";
  return aPath;
}

std::array<double, 6> bounds(const TopoDS_Shape& theShape)
{
  Bnd_Box aBox;
  BRepBndLib::Add(theShape, aBox);
  std::array<double, 6> aValues;
  aBox.Get(aValues[0], aValues[1], aValues[2], aValues[3], aValues[4], aValues[5]);
  return aValues;
}

double maxEdgeTolerance(const TopoDS_Shape& theShape)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anEdges;
  TopExp::MapShapes(theShape, TopAbs_EDGE, anEdges);
  double aMaximum = 0.0;
  for (int anIndex = 1; anIndex <= anEdges.Extent(); ++anIndex)
  {
    aMaximum = std::max(aMaximum, BRep_Tool::Tolerance(TopoDS::Edge(anEdges(anIndex))));
  }
  return aMaximum;
}
class Issue1761YFilletTest : public ::testing::TestWithParam<bool>
{
};
} // namespace

TEST_P(Issue1761YFilletTest, Face43And46RadiusThreeProducesValidSolid)
{
  TopoDS_Shape aInput;
  BRep_Builder aBuilder;
  ASSERT_TRUE(BRepTools::Read(aInput, fixturePath().ToCString(), aBuilder));
  ASSERT_TRUE(BRepCheck_Analyzer(aInput).IsValid());

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aFaces;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aAllEdges;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aSelected;
  TopExp::MapShapes(aInput, TopAbs_FACE, aFaces);
  TopExp::MapShapes(aInput, TopAbs_EDGE, aAllEdges);
  ASSERT_EQ(aFaces.Extent(), 133);
  ASSERT_EQ(aAllEdges.Extent(), 313);
  for (const int aFaceIndex : {43, 46})
  {
    for (TopExp_Explorer anExplorer(aFaces(aFaceIndex), TopAbs_EDGE); anExplorer.More();
         anExplorer.Next())
    {
      aSelected.Add(anExplorer.Current());
    }
  }
  ASSERT_EQ(aSelected.Extent(), 8);

  // Freeze the original issue's selected physical edges, not just the face numbers.
  const std::array<int, 8> anExpectedGlobalIndices = {44, 148, 149, 150, 151, 155, 47, 156};
  for (int anIndex = 1; anIndex <= aSelected.Extent(); ++anIndex)
  {
    EXPECT_EQ(aAllEdges.FindIndex(aSelected(anIndex)),
              anExpectedGlobalIndices[static_cast<size_t>(anIndex - 1)]);
  }

  BRepFilletAPI_MakeFillet aFillet(aInput);
  for (int anIteration = 0; anIteration < aSelected.Extent(); ++anIteration)
  {
    const int anIndex = GetParam() ? aSelected.Extent() - anIteration : anIteration + 1;
    aFillet.Add(3.0, TopoDS::Edge(aSelected(anIndex)));
  }
  aFillet.Build();
  ASSERT_TRUE(aFillet.IsDone());
  const TopoDS_Shape& aResult = aFillet.Shape();
  ASSERT_FALSE(aResult.IsNull());

  // IsDone alone was already true on the failing 4db baseline.
  ASSERT_TRUE(BRepCheck_Analyzer(aResult).IsValid());
  BRepAlgoAPI_Check anAdvancedCheck(aResult, true, true);
  ASSERT_TRUE(anAdvancedCheck.IsValid());

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aSolids;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(aResult, TopAbs_SOLID, aSolids);
  TopExp::MapShapes(aResult, TopAbs_EDGE, aResultEdges);
  ASSERT_EQ(aSolids.Extent(), 1);

  GProp_GProps anInputMass;
  GProp_GProps aResultMass;
  BRepGProp::VolumeProperties(aInput, anInputMass);
  BRepGProp::VolumeProperties(aResult, aResultMass);
  EXPECT_NEAR(anInputMass.Mass(), 2546681.4803456254, 0.01);
  EXPECT_NEAR(aResultMass.Mass(), 2546593.13937, 0.01);

  const auto anInputBounds = bounds(aInput);
  const auto aResultBounds = bounds(aResult);
  for (size_t anIndex = 0; anIndex < anInputBounds.size(); ++anIndex)
  {
    EXPECT_NEAR(aResultBounds[anIndex], anInputBounds[anIndex], 1.0e-5);
  }

  // The source already has 0.193296876 mm maximum edge tolerance.
  EXPECT_LE(maxEdgeTolerance(aResult), maxEdgeTolerance(aInput) + 1.0e-8);

  for (int anIndex = 1; anIndex <= aSelected.Extent(); ++anIndex)
  {
    const TopoDS_Edge aSeed = TopoDS::Edge(aSelected(anIndex));
    EXPECT_FALSE(aResultEdges.Contains(aSeed));
    bool hasRadiusThree = false;
    for (const TopoDS_Shape& aGenerated : aFillet.Generated(aSeed))
    {
      if (aGenerated.ShapeType() != TopAbs_FACE)
      {
        continue;
      }
      BRepAdaptor_Surface aSurface(TopoDS::Face(aGenerated));
      hasRadiusThree =
        hasRadiusThree
        || (aSurface.GetType() == GeomAbs_Cylinder
            && std::abs(aSurface.Cylinder().Radius() - 3.0) < 1.0e-6);
    }
    EXPECT_TRUE(hasRadiusThree) << "selected edge " << anIndex;
  }
}

INSTANTIATE_TEST_SUITE_P(ForwardAndReverse, Issue1761YFilletTest, ::testing::Bool());
