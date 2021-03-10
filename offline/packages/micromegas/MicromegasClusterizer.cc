/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasClusterizer.h"
#include "MicromegasDefs.h"
#include "CylinderGeomMicromegas.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>           // for PHG4CylinderGeom

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>                  // for TrkrCluster
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                     // for PHIODataNode
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>                   // for PHNodeIterator
#include <phool/PHObject.h>                         // for PHObject

#include <Eigen/Dense>

#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <cstdint>                                 // for uint16_t
#include <iterator>                                 // for distance
#include <map>                                      // for _Rb_tree_const_it...
#include <utility>                                  // for pair, make_pair
#include <vector>


namespace
{
  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

//_______________________________________________________________________________
MicromegasClusterizer::MicromegasClusterizer(const std::string &name, const std::string& detector)
  : SubsysReco(name)
  , m_detector( detector )
{}

//_______________________________________________________________________________
int MicromegasClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  assert( dstNode );

  // Create the Cluster node if missing
  auto trkrClusterContainer = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrClusterContainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto trkrNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if(!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    trkrClusterContainer = new TrkrClusterContainer();
    auto TrkrClusterContainerNode = new PHIODataNode<PHObject>(trkrClusterContainer, "TRKR_CLUSTER", "PHObject");
    trkrNode->addNode(TrkrClusterContainerNode);
  }

  // create cluster to hit association node, if missing
  auto trkrClusterHitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
  if(!trkrClusterHitAssoc)
  {
    PHNodeIterator dstiter(dstNode);
    auto trkrNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if(!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    trkrClusterHitAssoc = new TrkrClusterHitAssoc();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(trkrClusterHitAssoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    trkrNode->addNode(newNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________________________________________________
int MicromegasClusterizer::process_event(PHCompositeNode *topNode)
{

  // geometry
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  assert(geonode);

  // hitset container
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert( trkrhitsetcontainer );

  // cluster container
  auto trkrClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert( trkrClusterContainer );

  // cluster-hit association
  auto trkrClusterHitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  assert( trkrClusterHitAssoc );

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,
							 "ActsTrackingGeometry");
  assert( m_tGeometry );

  // loop over micromegas hitsets
  const auto hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId);
  for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
  {

    // get hitset, key and layer
    TrkrHitSet* hitset = hitset_it->second;
    const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;
    const auto layer = TrkrDefs::getLayer(hitsetkey);
    const auto tileid = MicromegasDefs::getTileId(hitsetkey);

    // get geometry object
    const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(geonode->GetLayerGeom(layer));
    assert(layergeom);

    /*
     * get segmentation type, layer thickness, strip length and pitch.
     * They are used to calculate cluster errors
     */
    const auto segmentation_type = layergeom->get_segmentation_type();
    const double thickness = layergeom->get_thickness();
    const double radius = layergeom->get_radius();
    const double pitch = layergeom->get_pitch();
    const double strip_length = layergeom->get_strip_length( tileid );

    // keep a list of ranges corresponding to each cluster
    using range_list_t = std::vector<TrkrHitSet::ConstRange>;
    range_list_t ranges;

    // loop over hits
    const auto hit_range = hitset->getHits();

    // keep track of first iterator of runing cluster
    auto begin = hit_range.first;

    // keep track of previous strip
    uint16_t previous_strip = 0;
    bool first = true;

    for( auto hit_it = hit_range.first; hit_it != hit_range.second; ++hit_it )
    {

      // get hit key
      const auto hitkey = hit_it->first;

      // get strip number
      const auto strip = MicromegasDefs::getStrip( hitkey );

      if( first )
      {

        previous_strip = strip;
        first = false;
        continue;

      } else if( strip - previous_strip > 1 ) {

        // store current cluster range
        ranges.push_back( std::make_pair( begin, hit_it ) );

        // reinitialize begin of next cluster range
        begin = hit_it;

      }

      // update previous strip
      previous_strip = strip;

    }

    // store last cluster
    if( begin != hit_range.second ) ranges.push_back( std::make_pair( begin, hit_range.second ) );

    // initialize cluster count
    int cluster_count = 0;

    // loop over found hit ranges and create clusters
    for( const auto& range : ranges )
    {

      // create cluster key and corresponding cluster
      const auto cluster_key = MicromegasDefs::genClusterKey( hitsetkey, cluster_count++ );
      auto cluster = (trkrClusterContainer->findOrAddCluster(cluster_key))->second;

      TVector3 world_coordinates;
      double weight_sum = 0;

      // needed for proper error calculation
      // it is either the sum over z, or phi, depending on segmentation
      double coord_sum = 0;
      double coordsquare_sum = 0;

      // loop over constituting hits
      for( auto hit_it = range.first; hit_it != range.second; ++hit_it )
      {
        // get hit key
        const auto hitkey = hit_it->first;
        const auto hit = hit_it->second;

        // associate cluster key to hit key
        trkrClusterHitAssoc->addAssoc(cluster_key, hitkey );

        // get strip number
        const auto strip = MicromegasDefs::getStrip( hitkey );

        // get adc, remove pedestal
        /* pedestal should be the same as the one used in PHG4MicromegasDigitizer */
        static constexpr double pedestal = 74.6;
        const double weight = double(hit->getAdc()) - pedestal;

        // get strip world coordinate and update relevant sums
        const auto strip_world_coordinate = layergeom->get_world_coordinate( tileid, strip );
        world_coordinates += strip_world_coordinate*weight;
        switch( segmentation_type )
        {
          case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
          {

            const auto rphi = radius*std::atan2( strip_world_coordinate.y(), strip_world_coordinate.x() );
            coord_sum += rphi*weight;
            coordsquare_sum += square(rphi)*weight;
            break;
          }

          case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
          {
            const auto z = strip_world_coordinate.z();
            coord_sum += z*weight;
            coordsquare_sum += square(z)*weight;
            break;
          }
        }

        weight_sum += weight;

      }

      // cluster position
      cluster->setPosition( 0, world_coordinates.x()/weight_sum );
      cluster->setPosition( 1, world_coordinates.y()/weight_sum );
      cluster->setPosition( 2, world_coordinates.z()/weight_sum );
      cluster->setGlobal();

      // dimension and error in r, rphi and z coordinates
      static const float invsqrt12 = 1./std::sqrt(12);
      static constexpr float error_scale_phi = 1.6;
      static constexpr float error_scale_z = 0.8;

      using matrix_t = Eigen::Matrix<float, 3, 3>;
      matrix_t dimension = matrix_t::Zero();
      matrix_t error = matrix_t::Zero();

      const auto size = std::distance( range.first, range.second );
      auto coord_cov = coordsquare_sum/weight_sum - square( coord_sum/weight_sum );
      auto coord_error_sq = coord_cov/weight_sum;
      switch( segmentation_type )
      {
        case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
        dimension(0,0) = square(0.5*thickness);
        dimension(1,1) = square(0.5*pitch*size);
        dimension(2,2) = square(0.5*strip_length);

        if( coord_error_sq == 0 ) coord_error_sq = square(pitch)/12;
        else coord_error_sq *= square(error_scale_phi);
        error(0,0) = square(thickness*invsqrt12);
        error(1,1) = coord_error_sq;
        error(2,2) = square(strip_length*invsqrt12);
        break;

        case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
        dimension(0,0) = square(0.5*thickness);
        dimension(1,1) = square(0.5*strip_length);
        dimension(2,2) = square(0.5*pitch*size);

        if( coord_error_sq == 0 ) coord_error_sq = square(pitch)/12;
        else coord_error_sq *= square(error_scale_z);
        error(0,0) = square(thickness*invsqrt12);
        error(1,1) = square(strip_length*invsqrt12);
        error(2,2) = coord_error_sq;
        break;
      }

      /// Add Acts and local information
           Acts::Vector3D globalPos(cluster->getX(), cluster->getY(), cluster->getZ());
          
      /// Get the surface key to find the surface from the map
      auto surface = getMmSurfaceFromCoords(topNode,
					    hitsetkey,
					    globalPos);
      
      Acts::Vector3D center = surface->center(m_tGeometry->geoContext)
	/ Acts::UnitConstants::cm;
      Acts::Vector3D normal = surface->normal(m_tGeometry->geoContext);
      
      double surfRadius = sqrt(center[0]*center[0] + center[1]*center[1]);
      double surfPhiCenter = atan2(center[1], center[0]);
      double surfRphiCenter = surfPhiCenter * surfRadius;
      double surfZCenter = center[2];
      
      double clusRadius = sqrt(cluster->getX() * cluster->getX() 
			       + cluster->getY() * cluster->getY());
      double clusphi = atan2(cluster->getY(), cluster->getX());
      double rClusPhi = clusRadius * clusphi;
      double zMm = globalPos(2);
      auto vecResult = surface->globalToLocal(m_tGeometry->geoContext, 
					      globalPos * Acts::UnitConstants::cm,
					      normal);
      Acts::Vector2D local2D;
      if(vecResult.ok())
	{
	  local2D = vecResult.value() / Acts::UnitConstants::cm;
	}
      else
	{
	  /// Otherwise use manual calculation, which is the same as Acts
	  local2D(0) = rClusPhi - surfRphiCenter;
	  local2D(1) = zMm - surfZCenter;
	}
  
      cluster->setLocalX(local2D(0));
      cluster->setLocalY(local2D(1));
      cluster->setActsSurface(surface);
      cluster->setActsLocalError(0,0, error(1,1));
      cluster->setActsLocalError(0,1, error(1,2));
      cluster->setActsLocalError(1,0, error(2,1));
      cluster->setActsLocalError(1,1,error(2,2));

      // rotate and save
      matrix_t rotation = matrix_t::Identity();
      const double phi = std::atan2(world_coordinates.y(), world_coordinates.x());
      const double cosphi = std::cos(phi);
      const double sinphi = std::sin(phi);
      rotation(0,0) = cosphi;
      rotation(0,1) = -sinphi;
      rotation(1,0) = sinphi;
      rotation(1,1) = cosphi;

      // rotate dimension and error
      dimension = rotation*dimension*rotation.transpose();
      error = rotation*error*rotation.transpose();

      // assign to cluster
      for( int i = 0; i<3; ++i )
        for( int j = 0; j<3; ++j )
      {
        cluster->setSize( i, j, dimension(i,j) );
        cluster->setError( i, j, error(i,j) );
      }

    }

  }

  // done
  return Fun4AllReturnCodes::EVENT_OK;
}

Surface MicromegasClusterizer::getMmSurfaceFromCoords(PHCompositeNode *topNode,
						      TrkrDefs::hitsetkey hitsetkey, 
						      Acts::Vector3D world)
{
  
  auto surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  if(!surfMaps)
    {
      std::cout << PHWHERE << "ActsSurfaceMaps not found on node tree!No TPOT clusters will be created."
		<< std::endl;
      return nullptr;
    }

  std::map<TrkrDefs::hitsetkey, std::vector<Surface>>::iterator mapIter;
  mapIter = surfMaps->mmSurfaceMap.find(hitsetkey);
  
  if(mapIter == surfMaps->mmSurfaceMap.end())
    {
      std::cout << PHWHERE 
		<< "Error: hitsetkey not found in clusterSurfaceMap, hitsetkey = "
		<< hitsetkey << std::endl;
      return nullptr;
    }

  double world_phi = atan2(world[1], world[0]);
  double world_z = world[2];
  
  std::vector<Surface> surf_vec = mapIter->second;
  unsigned int surf_index = 999;

  /// Get some geometry values from the geom builder for parsing surfaces
  double surfStepPhi = m_tGeometry->mmSurfStepPhi;
  double surfStepZ = m_tGeometry->mmSurfStepZ;

  for(unsigned int i=0;i<surf_vec.size(); ++i)
    {
      Surface this_surf = surf_vec[i];
  
      auto vec3d = this_surf->center(m_tGeometry->geoContext);
      std::vector<double> surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
      double surf_phi = atan2(surf_center[1], surf_center[0]);
      double surf_z = surf_center[2];

      /// Check if the cluster is geometrically within the surface boundaries
      /// The MMs surfaces span the entire length in z, so we don't divide
      /// by 2 in the z direction since the center of the surface is z=0
      bool withinPhi = world_phi >= surf_phi - surfStepPhi / 2.0
	&& world_phi < surf_phi + surfStepPhi / 2.0;
      bool withinZ = world_z > surf_z - surfStepZ 
	&& world_z < surf_z + surfStepZ;
      if( withinPhi && withinZ )
	{
	  surf_index = i;	  
	  break;
	}
    }
  if(surf_index == 999)
    {
      std::cout << PHWHERE 
		<< "Error: Micromegas surface index not defined, skipping cluster!" 
		<< std::endl;
    
      return nullptr;
    }
 
  return surf_vec[surf_index];

}
