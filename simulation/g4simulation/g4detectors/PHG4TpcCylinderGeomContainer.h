// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4TPCCYLINDERGEOMCONTAINER_H
#define G4DETECTORS_PHG4TPCCYLINDERGEOMCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for make_pair, pair

class PHG4TpcCylinderGeom;

class PHG4TpcCylinderGeomContainer: public PHObject
{
 public:
  using Map = std::map<int,PHG4TpcCylinderGeom *>;
  using Iterator = Map::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  PHG4TpcCylinderGeomContainer();
  ~PHG4TpcCylinderGeomContainer() override;

// from PHObject
  void identify(std::ostream& os = std::cout) const override;

  int AddLayerGeom(const int i, PHG4TpcCylinderGeom *mygeom);
  int AddLayerGeom(PHG4TpcCylinderGeom *mygeom);
  PHG4TpcCylinderGeom *GetLayerGeom(const int i);
  PHG4TpcCylinderGeom *GetFirstLayerGeom();
  int get_NLayers() const {return layergeoms.size();}
  ConstRange get_begin_end() const {return std::make_pair(layergeoms.begin(), layergeoms.end());}

 protected:
  Map layergeoms;
  float magfield;
  ClassDefOverride(PHG4TpcCylinderGeomContainer,1)
};

#endif
