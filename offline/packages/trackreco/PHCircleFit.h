#ifndef TRACKRECO_PHCIRCLEFIT_H
#define TRACKRECO_PHCIRCLEFIT_H


#include <fun4all/SubsysReco.h>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Units.hpp>

class SvtxTrackMap;
class TrkrClusterContainer;
class TrkrCluster;
class SvtxTrack;
class SvtxVertexMap;

class PHCircleFit : public SubsysReco
{

 public:
  
  PHCircleFit(const std::string& name = "PHCircleFit");

  ~PHCircleFit() {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode *topNode) override;

  void minLayer(int layer)
    { m_minLayer = layer; }
  void maxLayer(int layer)
    { m_maxLayer = layer; }
 
 private:

  int getNodes(PHCompositeNode *topNode);
  void circleFitTrack(SvtxTrack *track,
		      Acts::Vector3D& position,
		      Acts::Vector3D& momentum,
		      int& charge);
  std::vector<TrkrCluster*> getClusters(SvtxTrack *track);
  void circleFitByTaubin(std::vector<TrkrCluster*>& clusters,
			 double& R, double& X0, double& Y0);
  void findRoot(double& R, double& X0, double& Y0,
		double& x, double& y, const int& vertexId);
  int getCharge(const std::vector<TrkrCluster*>& clusters,
		const double& circPhi);
  void lineFit(const std::vector<TrkrCluster*>& clusters,
	       double& A, double& B);
  double evaluateZVertex(const int vtxid, double& m, double& B, 
			 double& x, double& y);
  double normPhi2Pi(const double& phi);


  SvtxTrackMap *m_trackMap = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;

  int m_minLayer = 0;
  int m_maxLayer = 53;

};


#endif //TRACKRECO_PHCIRCLEFIT_H
