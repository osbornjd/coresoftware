#include "PHActsInitialVertexFinder.h"
#include "ActsTransformations.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>
#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include <Acts/Vertexing/IterativeVertexFinder.hpp>
#include <Acts/Vertexing/LinearizedTrack.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexFinderConcept.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>
#include <Acts/Vertexing/ZScanVertexFinder.hpp>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <memory>

PHActsInitialVertexFinder::PHActsInitialVertexFinder(const std::string& name)
  : PHInitVertexing(name)
{}

int PHActsInitialVertexFinder::Setup(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  std::ostringstream name;
  name.str("");
  name << Name() <<"IVF.root";
  m_file = new TFile(name.str().c_str(),"RECREATE");
  dumhist = new TH1F("dumhist","dumhist",100,0,10);
  tree = new TTree("ivftree","tree with vertices");
  tree->Branch("m_vx", &m_vx,"m_vx/F");
  tree->Branch("m_vy", &m_vy,"m_vy/F");
  tree->Branch("m_vz", &m_vz,"m_vz/F");  
  tree->Branch("m_ntrk", &m_ntrk,"m_ntrk/F");
  tree->Branch("m_chi2",&m_chi2,"m_chi2/F");
  tree->Branch("m_ndf", &m_ndf,"m_ndf/F");
  tree->Branch("m_x", &m_x);
  tree->Branch("m_y", &m_y);
  tree->Branch("m_z", &m_z);  
  tree->Branch("m_px", &m_px);
  tree->Branch("m_py", &m_py);
  tree->Branch("m_pz", &m_pz);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::Process(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    std::cout << "PHActsInitialVertexFinder processing event " 
	      << m_event << std::endl;

  InitKeyMap keyMap;
  auto trackPointers = getTrackPointers(keyMap);

  auto vertices = findVertices(trackPointers);

  fillVertexMap(vertices, keyMap);

  /// Need to check that silicon stubs which may have been
  /// skipped over still have a vertex association
  checkTrackVertexAssociation();

  for(auto track : trackPointers)
    {
      delete track;
    }

  if(Verbosity() > 0)
    std::cout << "PHActsInitialVertexFinder processed event "
	      << m_event << std::endl;

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::End(PHCompositeNode *topNode)
{
  m_file->cd();
  dumhist->Write();
  tree->Write();
  m_file->Write();
  m_file->Close();
  std::cout << "Acts IVF used " << m_propagateTotals << " prop results out of "
	    << m_trackTotals << " total tracks " << std::endl;
  std::cout << "Acts IVF succeeded " << m_successFits 
	    << " out of " << m_totVertexFits << " total fits"
	    << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsInitialVertexFinder::checkTrackVertexAssociation()
{
  
  for(auto& [trackKey, track] : *m_trackMap)
    {
      /// If the track wasn't already given a vertex ID, it wasn't 
      /// included in the initial vertex finding due to Acts not liking
      /// tracks with large transverse position. So find the closest
      /// z vertex to it and assign it
      if(track->get_vertex_id() != UINT_MAX)
	continue;

      const auto trackZ = track->get_z();
      
      double closestVertZ = 9999;
      int vertId = -1;
      for(auto& [vertexKey, vertex] : *m_vertexMap)
	{
	  double dz = fabs(trackZ - vertex->get_z());

	  if(dz < closestVertZ) 
	    {
	      vertId = vertexKey;
	      closestVertZ = dz;
	    }

	}
      track->set_vertex_id(vertId);	

    }

}
void PHActsInitialVertexFinder::fillVertexMap(VertexVector& vertices,
					      InitKeyMap& keyMap)
{
  unsigned int vertexId = 0;
  if(vertices.size() != 0)
    m_vertexMap->clear();

  /// Create a fail safe for (e.g.) single particle events which 
  /// don't return a vertex
  if(vertices.size() == 0)
    {
      createDummyVertex();
      if(Verbosity() > 1)
	std::cout << "No vertices found. Adding a dummy vertex"
		  << std::endl;
      return;
    }

  for(auto vertex : vertices)
    {
      const auto &[chi2, ndf] = vertex.fitQuality();
      const auto numTracks = vertex.tracks().size();
      
      if(Verbosity() > 1)
	{
	  std::cout << "Found vertex at (" << vertex.position().x()
		    << ", " << vertex.position().y() << ", " 
		    << vertex.position().z() << ")" << std::endl;
	  std::cout << "Vertex has ntracks = " << numTracks
		    << " with chi2/ndf " << chi2 / ndf << std::endl;
	}

      /// Fill SvtxVertexMap
      #if __cplusplus < 201402L
      auto svtxVertex = boost::make_unique<SvtxVertex_v1>();
      #else
      auto svtxVertex = std::make_unique<SvtxVertex_v1>();
      #endif

      svtxVertex->set_x(vertex.position().x() / Acts::UnitConstants::cm);  
      svtxVertex->set_y(vertex.position().y() / Acts::UnitConstants::cm);
      svtxVertex->set_z(vertex.position().z() / Acts::UnitConstants::cm);
      for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	  svtxVertex->set_error(i, j,
				vertex.covariance()(i,j) 
				/ Acts::UnitConstants::cm2); 
	}
      }
      svtxVertex->set_chisq(chi2);
      svtxVertex->set_ndof(ndf);
      svtxVertex->set_t0(vertex.time());
      svtxVertex->set_id(vertexId);
          
      m_vx = svtxVertex->get_x();
      m_vy = svtxVertex->get_y();
      m_vz = svtxVertex->get_z();
      m_ntrk = numTracks;
      m_chi2 = chi2;
      m_ndf = ndf;
      m_x.clear();
      m_y.clear();
      m_z.clear();
      m_px.clear();
      m_py.clear();
      m_pz.clear();

      for(const auto track : vertex.tracks())
	{
	  const auto originalParams = track.originalParams;

	  const auto trackKey = keyMap.find(originalParams)->second;
	  svtxVertex->insert_track(trackKey);

	  /// Give the track the appropriate vertex id
	  const auto svtxTrack = m_trackMap->find(trackKey)->second;
	  m_x.push_back(svtxTrack->get_x());
	  m_y.push_back(svtxTrack->get_y());
	  m_z.push_back(svtxTrack->get_z());
	  m_px.push_back(svtxTrack->get_px());
	  m_py.push_back(svtxTrack->get_py());
	  m_pz.push_back(svtxTrack->get_pz());
	  if(Verbosity() > 3)
	    {   
	      svtxTrack->identify();
	      std::cout << "Updating track key " << trackKey << " with vertex "
			<< vertexId << std::endl;
	    }

	  svtxTrack->set_vertex_id(vertexId);
	}
      tree->Fill();
      m_vertexMap->insert(svtxVertex.release());

      ++vertexId;
    }
      
  return;
}

void PHActsInitialVertexFinder::createDummyVertex()
{

  /// If the Acts IVF finds 0 vertices, there weren't enough tracks
  /// for it to properly identify a good vertex. So just create
  /// a dummy vertex with large covariance for rest of track
  /// reconstruction to avoid seg faults

  #if __cplusplus < 201402L
  auto svtxVertex = boost::make_unique<SvtxVertex_v1>();
  #else
  auto svtxVertex = std::make_unique<SvtxVertex_v1>();
  #endif

  svtxVertex->set_x(0);  
  svtxVertex->set_y(0);
  svtxVertex->set_z(0);
  
  for(int i = 0; i < 3; ++i) 
    for(int j = 0; j < 3; ++j)
      {
	if( i == j)
	  svtxVertex->set_error(i, j, 10.); 
	else 
	  svtxVertex->set_error(i,j, 0);
      }
  float nan = NAN;
  svtxVertex->set_chisq(nan);
  svtxVertex->set_ndof(nan);
  svtxVertex->set_t0(nan);
  svtxVertex->set_id(0);

  m_vertexMap->insert(svtxVertex.release());
  
  for(auto& [key, track] : *m_trackMap)
    track->set_vertex_id(0);

}
 
VertexVector PHActsInitialVertexFinder::findVertices(TrackParamVec& tracks)
{

  m_totVertexFits++;
  
  /// Determine the input mag field type from the initial geometry
  /// and run the vertex finding with the determined mag field
  return std::visit([tracks, this](auto &inputField) {
      /// Setup aliases
      using InputMagneticField = 
	typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField = Acts::SharedBField<InputMagneticField>;
      
      using Stepper = Acts::EigenStepper<MagneticField>;
      using Propagator = Acts::Propagator<Stepper>;
      using TrackParameters = Acts::BoundTrackParameters;
      using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
      using VertexFitter = 
	Acts::FullBilloirVertexFitter<TrackParameters,Linearizer>;
      using ImpactPointEstimator = 
	Acts::ImpactPointEstimator<TrackParameters, Propagator>;
      using VertexSeeder = Acts::ZScanVertexFinder<VertexFitter>;
      using VertexFinder = 
	Acts::IterativeVertexFinder<VertexFitter, VertexSeeder>;
      using VertexFinderOptions = Acts::VertexingOptions<TrackParameters>;

      static_assert(Acts::VertexFinderConcept<VertexSeeder>,
		    "VertexSeeder does not fulfill vertex finder concept.");
      static_assert(Acts::VertexFinderConcept<VertexFinder>,
		    "VertexFinder does not fulfill vertex finder concept.");

      auto logLevel = Acts::Logging::FATAL;
      if(Verbosity() > 4)
	logLevel = Acts::Logging::VERBOSE;
      auto logger = Acts::getDefaultLogger("PHActsInitialVertexFinder", 
					   logLevel);
    
      MagneticField bField(inputField);
      auto propagator = std::make_shared<Propagator>(Stepper(bField));
      
      /// Setup vertex finder now
      typename VertexFitter::Config vertexFitterConfig;

      /// Vertex fitter seems to have no performance difference when
      /// iterating once vs. default of 5 times. Additionally, iterating
      /// more than once causes vertices with low numbers of tracks to
      /// fail fitting, causing an error to be thrown and 0 vertices 
      /// returned
      vertexFitterConfig.maxIterations = m_maxIterations;

      VertexFitter vertexFitter(std::move(vertexFitterConfig));
      
      typename Linearizer::Config linearizerConfig(bField, propagator);
      Linearizer linearizer(std::move(linearizerConfig));
      
      typename ImpactPointEstimator::Config ipEstConfig(bField, propagator);
      ImpactPointEstimator ipEst(std::move(ipEstConfig));
      
      typename VertexSeeder::Config seederConfig(ipEst);

      /// Don't weight track contribution by pT, since the momentum
      /// resolution of the silicon seeds is poor
      if(m_disableWeights)
	seederConfig.disableAllWeights = true;
      VertexSeeder seeder(std::move(seederConfig));
      
      typename VertexFinder::Config finderConfig(std::move(vertexFitter), 
						 std::move(linearizer),
						 std::move(seeder), ipEst);
      finderConfig.maxVertices = m_maxVertices;
      finderConfig.reassignTracksAfterFirstFit = true;
      VertexFinder finder(finderConfig, std::move(logger));

      typename VertexFinder::State state(m_tGeometry->magFieldContext);
      VertexFinderOptions finderOptions(m_tGeometry->geoContext,
					m_tGeometry->magFieldContext);
  
      auto result = finder.find(tracks, finderOptions, state);
    
      VertexVector vertexVector;

      if(result.ok())
	{
	  m_successFits++;

	  auto vertexCollection = *result;
	  
	  if(Verbosity() > 1)
	    {
	      std::cout << "Acts IVF found " << vertexCollection.size()
			<< " vertices in event" << std::endl;
	    }
	  
	  for(const auto& vertex : vertexCollection) 
	    {
	      vertexVector.push_back(vertex);
	    }
	}
      else
	{
	  if(Verbosity() > 0)
	    {
	      std::cout << "Acts initial vertex finder returned error: " 
			<< result.error().message() << std::endl;
	      std::cout << "Track positions IVF used are : " << std::endl;
	      for(const auto track : tracks)
		{
		  const auto position = track->position(m_tGeometry->geoContext);
		  std::cout << "(" << position(0) << ", " << position(1)
			    << ", " << position(2) << std::endl;
		}
	    }
	}
    
      return vertexVector;
      
    } /// end lambda
    , m_tGeometry->magField
    ); /// end std::visit call

}
TrackParamVec PHActsInitialVertexFinder::getTrackPointers(InitKeyMap& keyMap)
{
  TrackParamVec tracks;
  
  for(auto& [key,track] : *m_trackMap)
    {
      if(Verbosity() > 3)
	{
	  std::cout << "Adding track seed to vertex finder " 
		    << std::endl;
	  track->identify();
	}

      const Acts::Vector4D stubVec(
                  track->get_x() * Acts::UnitConstants::cm,
		  track->get_y() * Acts::UnitConstants::cm,
		  track->get_z() * Acts::UnitConstants::cm,
		  10 * Acts::UnitConstants::ns);
     
      const Acts::Vector3D stubMom(track->get_px(),
				   track->get_py(),
				   track->get_pz());
      const int trackQ = track->get_charge() * Acts::UnitConstants::e;
      const double p = track->get_p();
      
      /// Make a dummy loose covariance matrix for Acts that corresponds
      /// to the silicon seed resolutions
      Acts::BoundSymMatrix cov;
      if(m_resetTrackCovariance && m_siliconSeeds)
	{
	  cov << 5000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
	         0., 900 * Acts::UnitConstants::um, 0., 0., 0., 0.,
	         0., 0., 0.005, 0., 0., 0.,
	         0., 0., 0., 0.001, 0., 0.,
	         0., 0., 0., 0., 0.3 , 0.,
	         0., 0., 0., 0., 0., 1.;
	}
      else if (m_resetTrackCovariance && !m_siliconSeeds)
	{
	  cov << 1000 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
	         0., 900 * Acts::UnitConstants::um, 0., 0., 0., 0.,
	         0., 0., 0.005, 0., 0., 0.,
	         0., 0., 0., 0.001, 0., 0.,
	         0., 0., 0., 0., 0.00005 , 0.,
	         0., 0., 0., 0., 0., 1.;
	}
      else 
	{
	  ActsTransformations transform;
	  transform.setVerbosity(Verbosity());
	  cov = transform.rotateSvtxTrackCovToActs(track,
						   m_tGeometry->geoContext);
	}
   
      const Acts::CurvilinearTrackParameters cParam(stubVec, stubMom,
						    trackQ / p, cov);
      auto result = propagateTrack(cParam);
      if(result.ok())
	{
	  std::unique_ptr<const Acts::BoundTrackParameters> params = std::move(*result);
	  /// The vertex finder requires raw pointers (not sure why)
	  /// so convert the unique_ptr to raw
	  const auto paramRawPtr = new Acts::BoundTrackParameters(*params);
	  tracks.push_back(paramRawPtr);
	  keyMap.insert(std::make_pair(paramRawPtr, key));
	  m_propagateTotals++;
	}
      else
	{

	  /// If the propagation doesn't work, we'll just add the
	  /// track on its own surface
	  /// Make a dummy perigeee surface to bind the track to
	  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
		 Acts::Vector3D(track->get_x() * Acts::UnitConstants::cm,
				track->get_y() * Acts::UnitConstants::cm,
				track->get_z() * Acts::UnitConstants::cm));
	  const auto paramPtr = new Acts::BoundTrackParameters(
				    perigee, m_tGeometry->geoContext,
				    stubVec, stubMom, p, trackQ, cov);
	  tracks.push_back(paramPtr);
	  keyMap.insert(std::make_pair(paramPtr, key));
	}
      m_trackTotals++;
    }

  return tracks;
}
BoundTrackParamPtrResult PHActsInitialVertexFinder::propagateTrack(
			 const Acts::CurvilinearTrackParameters param)
{
  /**
   * This function propagates the determined track parameters to a
   * perigee surface centered at 0,0,0 so that all parameters are
   * defined in the same way wrt to the same surface before going 
   * into the initial vertex finder. Nominally we would do this wrt
   * the initial vertex estimate, but of course we don't know that yet.
   */
  
  return std::visit([param, this] (auto &&inputField) -> BoundTrackParamPtrResult {

      using InputMagneticField = typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField = Acts::SharedBField<InputMagneticField>;
      using Stepper = Acts::EigenStepper<MagneticField>;
      using Propagator = Acts::Propagator<Stepper>;
      
      
      MagneticField field(inputField);
      Stepper stepper(field);
      Propagator propagator(stepper);
      
      Acts::Logging::Level logLevel = Acts::Logging::FATAL;
      if(Verbosity() > 3)
	logLevel = Acts::Logging::VERBOSE;
      
      auto logger = Acts::getDefaultLogger("PHActsInitialVertexFinder",
					   logLevel);
      Acts::PropagatorOptions<> options(m_tGeometry->geoContext,
					m_tGeometry->magFieldContext,
					Acts::LoggerWrapper{*logger});
      auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
			            Acts::Vector3D(0.,0.,0.));

      auto result = propagator.propagate(param, *perigee,
					 options);

      if(result.ok())
	{
	  return std::move((*result).endParameters);
	}
      else 
	{
	std::cout << "IVF propagation error... couldn't transport parameters to perigee" << std::endl;
	return result.error();
      }
      
    }, m_tGeometry->magField);


}
int PHActsInitialVertexFinder::getNodes(PHCompositeNode *topNode)
{

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_svtxTrackMapName.c_str());
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "No " << m_svtxTrackMapName.c_str() 
		<< " on node tree, bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

   m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, 
							 "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::createNodes(PHCompositeNode *topNode)
{
  
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  /// Check that it is there
  if (!dstNode)
  {
    std::cerr << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in PHActsInitialVertexFinder::createNodes");
  }

  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode,
						  m_svtxVertexMapName.c_str());
  
  if(!m_vertexMap)
    {
      m_vertexMap = new SvtxVertexMap_v1;
      PHIODataNode<PHObject>* vertexNode = new PHIODataNode<PHObject>( 
		   m_vertexMap, m_svtxVertexMapName.c_str(),"PHObject");

      svtxNode->addNode(vertexNode);

    }


  return Fun4AllReturnCodes::EVENT_OK;
}
