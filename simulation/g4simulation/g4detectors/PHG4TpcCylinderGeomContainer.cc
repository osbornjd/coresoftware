#include "PHG4TpcCylinderGeomContainer.h"
#include "PHG4TpcCylinderGeom.h"
#include <cmath>

using namespace std;

PHG4TpcCylinderGeomContainer::PHG4TpcCylinderGeomContainer():
  magfield(NAN)
{}

PHG4TpcCylinderGeomContainer::~PHG4TpcCylinderGeomContainer()
{
  while(layergeoms.begin() != layergeoms.end())
    {
      delete layergeoms.begin()->second;
      layergeoms.erase(layergeoms.begin());
    }
  return;
}

void
PHG4TpcCylinderGeomContainer::identify(std::ostream& os) const
{
  os << "mag field: " << magfield << endl;
  os << "number of layers: " << layergeoms.size() << endl;
  map<int,PHG4TpcCylinderGeom *>::const_iterator iter;
  for (iter=layergeoms.begin(); iter != layergeoms.end(); ++iter)
    {
      (iter->second)->identify(os);
    }

  return;
}

int
PHG4TpcCylinderGeomContainer::AddLayerGeom(const int i, PHG4TpcCylinderGeom *mygeom)
{
  if (layergeoms.find(i) != layergeoms.end())
    {
      cout << "layer " << i << " already added to PHCylinderGeomContainer" << endl;
      return -1;
    }
  mygeom->set_layer(i);
  layergeoms[i] = mygeom;
  return 0;
}

int
PHG4TpcCylinderGeomContainer::AddLayerGeom(PHG4TpcCylinderGeom *mygeom)
{
  int layer = mygeom->get_layer();
  if (layergeoms.find(layer) != layergeoms.end())
    {
      cout << "layer " << layer << " already added to PHCylinderGeomContainer" << endl;
      return -1;
    }
  layergeoms[layer] = mygeom;
  return 0;
}

PHG4TpcCylinderGeom *
PHG4TpcCylinderGeomContainer::GetLayerGeom(const int i)
{
  map<int,PHG4TpcCylinderGeom *>::const_iterator iter = layergeoms.find(i);
  if (iter !=  layergeoms.end())
    {
      return iter->second;
    }
  cout << "Could not locate layer " << i << " in PHG4TpcCylinderGeomContainer" << endl;
  return nullptr;
}

PHG4TpcCylinderGeom *
PHG4TpcCylinderGeomContainer::GetFirstLayerGeom()
{
  if (layergeoms.empty())
    {
      return nullptr;
    }
  return layergeoms.begin()->second;
}

