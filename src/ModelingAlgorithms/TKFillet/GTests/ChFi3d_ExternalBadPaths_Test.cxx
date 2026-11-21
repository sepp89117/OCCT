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
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_Shape.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>

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

  BRepFilletAPI_MakeFillet fillet(input, ChFi3d_Rational);
  fillet.SetParams(0.01, 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-5, 1.0e-3);
  fillet.SetContinuity(GeomAbs_C1, 0.01);
  fillet.Add(2.2, seed); // The official script's second, formerly BAD construction.
  EXPECT_NO_THROW(fillet.Build());
  ASSERT_TRUE(fillet.IsDone());
  ASSERT_GT(fillet.NbContours(), 0);
  for (int contour = 1; contour <= fillet.NbContours(); ++contour)
    EXPECT_NEAR(fillet.Radius(contour), 2.2, 1.0e-9);

  const TopoDS_Shape& result = fillet.Shape();
  ASSERT_FALSE(result.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(result).IsValid());
  BRepAlgoAPI_Check booleanCheck(result, true, true);
  ASSERT_TRUE(booleanCheck.IsValid());
  int radiusFaceCount = 0;
  for (TopExp_Explorer it(result, TopAbs_FACE); it.More(); it.Next())
  {
    const BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
    if (surface.GetType() == GeomAbs_Cylinder
        && std::abs(surface.Cylinder().Radius() - 2.2) < 1.0e-6)
      ++radiusFaceCount;
  }
  EXPECT_EQ(radiusFaceCount, 1);
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





namespace
{
void checkOfficialBadPath(const char* fixture,
                          double radius,
                          std::initializer_list<int> edgeIndices,
                          double expectedZMax,
                          double expectedTorusMinorRadius = 0.0,
                          std::initializer_list<int> historySeamIndices = {})
{
  TCollection_AsciiString path(__FILE__);
  const int slash = std::max(path.SearchFromEnd("/"), path.SearchFromEnd("\\"));
  path.Trunc(slash);
  path += "data/";
  path += fixture;
  TopoDS_Shape input;
  BRep_Builder builder;
  ASSERT_TRUE(BRepTools::Read(input, path.ToCString(), builder));
  ASSERT_TRUE(BRepCheck_Analyzer(input).IsValid());
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
  TopExp::MapShapes(input, TopAbs_EDGE, edges);

  BRepFilletAPI_MakeFillet fillet(input, ChFi3d_Rational);
  fillet.SetParams(0.01, 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-5, 1.0e-3);
  fillet.SetContinuity(GeomAbs_C1, 0.01);
  for (int index : edgeIndices)
  {
    ASSERT_GE(index, 1);
    ASSERT_LE(index, edges.Extent());
    fillet.Add(radius, TopoDS::Edge(edges(index)));
  }
  EXPECT_NO_THROW(fillet.Build());
  ASSERT_TRUE(fillet.IsDone());
  EXPECT_EQ(fillet.NbContours(), static_cast<int>(edgeIndices.size()));
  for (int contour = 1; contour <= fillet.NbContours(); ++contour)
    EXPECT_NEAR(fillet.Radius(contour), radius, 1.0e-9);
  for (int index : edgeIndices)
  {
    const TopoDS_Edge selected = TopoDS::Edge(edges(index));
    EXPECT_TRUE(fillet.IsDeleted(selected));
    int generatedFaces = 0;
    for (const TopoDS_Shape& generated : fillet.Generated(selected))
    {
      if (generated.ShapeType() != TopAbs_FACE)
        continue;
      ++generatedFaces;
      const BRepAdaptor_Surface surface(TopoDS::Face(generated));
      if (surface.GetType() == GeomAbs_Cylinder)
        EXPECT_NEAR(surface.Cylinder().Radius(), radius, 1.0e-6);
      else if (surface.GetType() == GeomAbs_Torus)
        EXPECT_NEAR(surface.Torus().MinorRadius(), radius, 1.0e-6);
      else
        ADD_FAILURE() << "Selected edge generated a non-analytic fillet surface";
    }
    EXPECT_EQ(generatedFaces, 1);
  }

  const TopoDS_Shape& result = fillet.Shape();
  ASSERT_FALSE(result.IsNull());
  ASSERT_TRUE(BRepCheck_Analyzer(result).IsValid());
  BRepAlgoAPI_Check booleanCheck(result, true, true);
  ASSERT_TRUE(booleanCheck.IsValid());
  // Every retained part of an original technical seam must be in the output.
  // A correct-looking solid with a stale split here breaks feature history.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> resultEdges;
  TopExp::MapShapes(result, TopAbs_EDGE, resultEdges);
  for (int index : historySeamIndices)
  {
    ASSERT_GE(index, 1);
    ASSERT_LE(index, edges.Extent());
    const TopoDS_Edge originalSeam = TopoDS::Edge(edges(index));
    const NCollection_List<TopoDS_Shape>& modified = fillet.Modified(originalSeam);
    EXPECT_EQ(modified.Extent(), 1);
    for (const TopoDS_Shape& part : modified)
    {
      EXPECT_EQ(part.ShapeType(), TopAbs_EDGE);
      EXPECT_TRUE(resultEdges.Contains(part));
    }
  }
  int solidCount = 0;
  for (TopExp_Explorer it(result, TopAbs_SOLID); it.More(); it.Next())
    ++solidCount;
  EXPECT_EQ(solidCount, 1);
  GProp_GProps originalVolume, resultVolume;
  BRepGProp::VolumeProperties(input, originalVolume);
  BRepGProp::VolumeProperties(result, resultVolume);
  EXPECT_GT(resultVolume.Mass(), 0.0);
  EXPECT_GT(std::abs(resultVolume.Mass() - originalVolume.Mass()), 1.0e-5);
  // The original OCC602 planar cap surfaces have loose analytical bounds.
  // Actual extremities are represented by the body's vertices.
  double zmin = std::numeric_limits<double>::infinity();
  double zmax = -std::numeric_limits<double>::infinity();
  for (TopExp_Explorer it(result, TopAbs_VERTEX); it.More(); it.Next())
  {
    const double z = BRep_Tool::Pnt(TopoDS::Vertex(it.Current())).Z();
    zmin = std::min(zmin, z);
    zmax = std::max(zmax, z);
  }
  EXPECT_NEAR(zmin, 0.0, 0.01);
  EXPECT_NEAR(zmax, expectedZMax, 0.01);

  if (expectedTorusMinorRadius > 0.0)
  {
    int torusCount = 0;
    for (TopExp_Explorer it(result, TopAbs_FACE); it.More(); it.Next())
    {
      const BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
      if (surface.GetType() != GeomAbs_Torus)
        continue;
      EXPECT_NEAR(surface.Torus().MinorRadius(), expectedTorusMinorRadius, 1.0e-6);
      ++torusCount;
    }
    EXPECT_EQ(torusCount, 1);

    auto capArea = [](const TopoDS_Shape& body, double zCap)
    {
      for (TopExp_Explorer it(body, TopAbs_FACE); it.More(); it.Next())
      {
        const TopoDS_Face face = TopoDS::Face(it.Current());
        const BRepAdaptor_Surface surface(face);
        if (surface.GetType() != GeomAbs_Plane
            || std::abs(surface.Plane().Location().Z() - zCap) > 1.0e-6)
          continue;
        GProp_GProps area;
        BRepGProp::SurfaceProperties(face, area);
        return area.Mass();
      }
      return -1.0;
    };
    for (const double zCap : {0.0, expectedZMax})
    {
      const double originalCap = capArea(input, zCap);
      const double resultCap = capArea(result, zCap);
      ASSERT_GT(originalCap, 0.0);
      EXPECT_NEAR(resultCap, originalCap, 1.0e-5);
    }
  }
}
}

TEST(ChFi3dExternalBadPaths, OCC602SeamFilletProducesValidSolid)
{
  // Exact official fixture SHA-256:
  // 693F7A5D8A0A78871402788CFFD04B9D2D7D8F0B9A6E8A3CDC6D603388123ED7.
  checkOfficialBadPath("OCC602.brep", 20.0, {9}, 220.0, 20.0);
}

TEST(ChFi3dExternalBadPaths, OCC26655TwelveEdgesProduceValidSolid)
{
  // Exact official fixture SHA-256:
  // C394355D83021B35FA874C397A028B55119BD7792BF22E93760B697ECFE1B8F0.
  checkOfficialBadPath("bug26655_OCC50.brep", 1.0,
                       {20, 10, 5, 17, 3, 4, 1, 2, 16, 11, 6, 7}, 30.0, 0.0,
                       {8, 12, 14, 18});
}
