#ifndef TRACKRECO_PHACTSINITIALVERTEXFINDER_H
#define TRACKRECO_PHACTSINITIALVERTEXFINDER_H

#include "PHInitVertexing.h"
#include <trackbase/ActsTrackingGeometry.h>

#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Result.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <TH1.h>
#include <TTree.h>
#include <TH2.h>
#include <TFile.h>
#include <sstream>
#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;

using VertexVector = std::vector<Acts::Vertex<Acts::BoundTrackParameters>>;

using BoundTrackParamPtr = 
  std::unique_ptr<const Acts::BoundTrackParameters>;
using TrackParamVec = std::vector<const Acts::BoundTrackParameters*>;

using InitKeyMap = std::map<const Acts::BoundTrackParameters*, const unsigned int>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;

class PHActsInitialVertexFinder: public PHInitVertexing
{
 public: 
  PHActsInitialVertexFinder(const std::string& name="PHActsInitialVertexFinder");
  virtual ~PHActsInitialVertexFinder() {}

  void setMaxVertices(const int maxVertices)
  { m_maxVertices = maxVertices;}

  void setSvtxVertexMapName(const std::string& name)
  { m_svtxVertexMapName = name; }
  
  void setSvtxTrackMapName(const std::string& name)
  { m_svtxTrackMapName = name; }

  void disablePtWeights(const bool weight)
  { m_disableWeights = weight; }

  void resetTrackCovariance(const bool initial)
  { m_resetTrackCovariance = initial; }
  
  void setSiliconSeeds(const bool silicon)
  { m_siliconSeeds = silicon; }
  
  void setMaxIterations(const int iter)
  { m_maxIterations = iter; }

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  
  TrackParamVec getTrackPointers(InitKeyMap& keyMap);
  VertexVector findVertices(TrackParamVec& tracks);
  void fillVertexMap(VertexVector& vertices, InitKeyMap& keyMap);
  void createDummyVertex();
  void checkTrackVertexAssociation();

  BoundTrackParamPtrResult propagateTrack(const Acts::CurvilinearTrackParameters param);
  
  int m_maxVertices = 5;
  int m_maxIterations = 1;
  int m_event = 0;
  unsigned int m_totVertexFits = 0;
  unsigned int m_successFits = 0;

  std::string m_svtxTrackMapName = "SvtxSiliconTrackMap";
  std::string m_svtxVertexMapName = "SvtxVertexMap";

  bool m_resetTrackCovariance = true;
  bool m_disableWeights = true;
  bool m_siliconSeeds = true;
  
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

  TFile *m_file= nullptr;
  TTree *tree = nullptr;
  float m_vx, m_vy, m_vz, m_ntrk, m_chi2, m_ndf;
  std::vector<float> m_x, m_y, m_z, m_px, m_py, m_pz;
  TH1* dumhist = nullptr;

  int m_trackTotals = 0;
  int m_propagateTotals = 0;
};


#endif
