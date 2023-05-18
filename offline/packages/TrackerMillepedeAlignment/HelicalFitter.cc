#include "HelicalFitter.h"

#include "Mille.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterContainer.h>   
#include <trackbase/TrkrClusterCrossingAssoc.h>   
#include <trackbase/alignmentTransformationContainer.h>

#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>

using namespace std;

//____________________________________________________________________________..
HelicalFitter::HelicalFitter(const std::string &name):
  SubsysReco(name)
  , PHParameterInterface(name)
  , _mille(nullptr)
{
  InitializeParameters();
}

//____________________________________________________________________________..
HelicalFitter::~HelicalFitter()
{

}

//____________________________________________________________________________..
int HelicalFitter::InitRun(PHCompositeNode *topNode)
{
  UpdateParametersWithMacro();

   int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output data file
  if(test_output)
    {
      _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  }
 else
   {
     _mille = new Mille(data_outfilename.c_str()); 
   }

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename);
  steering_file << data_outfilename << std::endl;
  steering_file.close();

  _file = new TFile(TString(Name() + "_hf.root"),"recreate");
  _residuallayerx = new TH2F("residuallayerx",";x residual [mm]; layer",
			     1000,-1,1,60,0,60);
  _residuallayerz = new TH2F("residuallayerz",";z residual [mm]; layer",
			     1000,-1,1,60,0,60);
  _residuallayerxpull = new TH2F("residuallayerxpull",";x residual pull; layer",
				 1000,-10, 10,60,0,60);
  _residuallayerzpull = new TH2F("residuallayerzpull",";z residual pull; layer",
				 1000,-10, 10,60,0,60);

  // print grouping setup to log file:
  std::cout << "MakeMilleFiles::InitRun: Surface groupings are silicon " << si_grp << " tpc " << tpc_grp << " mms " << mms_grp << std::endl; 

  return ret;
}

//_____________________________________________________________________
void HelicalFitter::SetDefaultParameters()
{


  return;
}

//____________________________________________________________________________..
int HelicalFitter::process_event(PHCompositeNode*)
{
  // _track_map_tpc contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds

  if(Verbosity() > 0)
    cout << PHWHERE 
	 << " TPC seed map size " << _track_map_tpc->size() 
	 << " Silicon seed map size "  << _track_map_silicon->size() 
	 << endl;

  if(_track_map_silicon->size() == 0 && _track_map_tpc->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;

  // Decide whether we want to make a helical fit for silicon or TPC
  unsigned int maxtracks = 0; 
  if(fittpc) { maxtracks =  _track_map_tpc->size();  }
  if(fitsilicon)  { maxtracks =  _track_map_silicon->size(); }
  for(unsigned int trackid = 0; trackid < maxtracks; ++trackid)
    {
      TrackSeed* tracklet = nullptr;
      if(fitsilicon) {  tracklet = _track_map_silicon->get(trackid); }
      else if(fittpc) {  tracklet = _track_map_tpc->get(trackid);	 }
      if(!tracklet) { trackid++;  continue; }

      std::vector<Acts::Vector3> global_vec;
      std::vector<TrkrDefs::cluskey> cluskey_vec;

      // Get a vector of cluster keys from the tracklet  
      getTrackletClusterList(tracklet, cluskey_vec);
      // store cluster global positions in a vector
      TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);   

      correctTpcGlobalPositions( global_vec, cluskey_vec);

      std::vector<float> fitpars =  TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit`

      if(fitpars.size() == 0) continue;  // discard this track, not enough clusters to fit

      if(Verbosity() > 1)  
	{ std::cout << " Track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
		 << " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }

      // if a full track is requested, get the silicon clusters too and refit
      if(fittpc && fitfulltrack)
	{
	  // this associates silicon clusters and adds them to the vectors
	  unsigned int nsilicon = TrackFitUtils::addSiliconClusters(fitpars, dca_cut, _tGeometry, _cluster_map, global_vec, cluskey_vec);
	  if(nsilicon < 3) continue;  // discard this TPC seed, did not get a good match to silicon

	  // fit the full track now
	  fitpars.clear();
	  fitpars =  TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit
	  if(fitpars.size() == 0) continue;  // discard this track, fit failed

	  if(Verbosity() > 1)  
	    { std::cout << " Full track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
					   << " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }
	} 

      // get the residuals and derivatives for all clusters

      for(unsigned int ivec=0;ivec<global_vec.size(); ++ivec)
	{

	  auto global = global_vec[ivec];
	  auto cluskey = cluskey_vec[ivec];
	  auto cluster = _cluster_map->findCluster(cluskey);
	  if(!cluster) { continue;}

	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);

	  // What we need now is to find the point on the surface at which the helix would intersect
	  // If we have that point, we can transform the fit back to local coords
	  // we have fitpars for the helix, and the cluster key - from which we get the surface

	  Surface surf = _tGeometry->maps().getSurface(cluskey, cluster);
	  Acts::Vector3 fitpoint = get_helix_surface_intersection(surf, fitpars, global);

	  // fitpoint is the point where the helical fit intersects the plane of the surface
	  // Now transform the helix fitpoint to local coordinates to compare with cluster local coordinates
	  Acts::Vector3 fitpoint_local = surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (fitpoint *  Acts::UnitConstants::cm);
	  fitpoint_local /= Acts::UnitConstants::cm;

	  auto xloc = cluster->getLocalX();  // in cm
	  auto zloc = cluster->getLocalY();	  
	  if(trkrid == TrkrDefs::tpcId) { zloc = convertTimeToZ(cluskey, cluster); }
	  //float yloc = 0.0;   // Because the fitpoint is on the surface, y will always be zero in local coordinates

	  Acts::Vector2 residual(xloc - fitpoint_local(0), zloc - fitpoint_local(1)); 

	  unsigned int layer = TrkrDefs::getLayer(cluskey_vec[ivec]);	  
	  float phi =  atan2(global(1), global(0));
	  float beta =  atan2(global(2), sqrt(pow(global(0),2) + pow(global(1),2)));
	  _residuallayerx->Fill(residual(0) * 10. , layer);
	  _residuallayerz->Fill(residual(1)*10., layer);
	  _residuallayerxpull->Fill(residual(0) / cluster->getRPhiError(),
				    layer);
	  _residuallayerzpull->Fill(residual(1) / cluster->getZError(), layer);

	  if(Verbosity() > 1) {
	  Acts::Vector3 loc_check =  surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (global *  Acts::UnitConstants::cm);
	  loc_check /= Acts::UnitConstants::cm;
	  std::cout << "    layer " << layer << std::endl
		    << " cluster global " << global(0) << " " << global(1) << " " << global(2) << std::endl
		    << " fitpoint " << fitpoint(0) << " " << fitpoint(1) << " " << fitpoint(2) << std::endl
		    << " fitpoint_local " << fitpoint_local(0) << " " << fitpoint_local(1) << " " << fitpoint_local(2) << std::endl  
		    << " cluster local x " << cluster->getLocalX() << " cluster local y " << cluster->getLocalY() << std::endl
		    << " cluster global to local x " << loc_check(0) << " local y " << loc_check(1) << "  local z " << loc_check(2) << std::endl
		    << " cluster local residual x " << residual(0) << " cluster local residual y " <<residual(1) << std::endl;
	  }

	  // need standard deviation of measurements
	  Acts::Vector2 clus_sigma = getClusterError(cluster, cluskey, global);
	  if(isnan(clus_sigma(0)) || isnan(clus_sigma(1)))  { continue; }

	  int glbl_label[AlignmentDefs::NGL];
	  if(layer < 7)
	    {
	      AlignmentDefs::getSiliconGlobalLabels(surf, glbl_label, si_grp);
	    }
	  else if (layer < 55)
	    {
	      AlignmentDefs::getTpcGlobalLabels(surf, glbl_label, tpc_grp);
	    }
	  else if (layer < 57)
	    {
	      AlignmentDefs::getMMGlobalLabels(surf, glbl_label, mms_grp);
	    }

	  // These derivatives are for the local parameters
	  // The angleDerivs dimensions are [alpha/beta/gamma](x/y/z)
	  //std::vector<Acts::Vector3> angleDerivs = getDerivativesAlignmentAngles(global, cluskey, cluster); 
	  //std::vector<Acts::Vector3> translDerivs = getDerivativesAlignmentTranslations(global, cluskey, cluster);

	  float lcl_derivativeX[AlignmentDefs::NLC];
	  float lcl_derivativeY[AlignmentDefs::NLC];

	  getLocalDerivativesXY(surf, global, fitpars, lcl_derivativeX, lcl_derivativeY, layer);

	  float glbl_derivativeX[AlignmentDefs::NGL];
	  float glbl_derivativeY[AlignmentDefs::NGL];
	  getGlobalDerivativesXY(surf, global, fitpoint, fitpars, glbl_derivativeX, glbl_derivativeY, layer);

	  for(unsigned int i = 0; i < AlignmentDefs::NGL; ++i) 
	    {
	      if( is_layer_param_fixed(layer, i) || is_layer_fixed(layer) )
		{
		  glbl_derivativeX[i] = 0;
		  glbl_derivativeY[i] = 0;
		}

	      if(trkrid == TrkrDefs::tpcId)
		{
		  unsigned int sector = TpcDefs::getSectorId(cluskey_vec[ivec]);	  
		  unsigned int side = TpcDefs::getSide(cluskey_vec[ivec]);	  
		  if(is_tpc_sector_fixed(layer, sector, side))
		    {
		      //if(i==0) std::cout << " param " << i << " layer " << layer << " sector " << sector << " side " << side << std::endl;
		      glbl_derivativeX[i] = 0;
		      glbl_derivativeY[i] = 0;
		    }
		}
	    }

	  // Add the measurement separately for each coordinate direction to Mille
	  // set the derivatives non-zero only for parameters we want to be optimized
	  // local parameter numbering is arbitrary:
	  float errinf = 1.0;

	  if(_layerMisalignment.find(layer) != _layerMisalignment.end())
	    {
	      errinf = _layerMisalignment.find(layer)->second;
	    }

	  // provides output that can be grep'ed to make plots of input to mille
	  if(Verbosity() > 1)
	    {
	      //if(layer < 7)
		{
		  // radius = fitpars[0],  X0 = fitpars[1],  Y0 = fitpars[2], zslope = fitpars[3], Z0  = fitpars[4] 
		  std::cout << "Local residualsX: layer " << layer << " phi " << phi * 180 / M_PI << " beta " << beta * 180.90 / M_PI
			    << " dxloc " << residual(0) << " error " << clus_sigma(0)  << " inflation_factor " << errinf
			    << " xloc " << xloc << " fitxloc " << fitpoint_local(0) 
			    << " zglob " << global(2) << " fitzglob " << fitpoint(2) 
			    << " xglob " << global(0) << " fitxglob " << fitpoint(0)
			    << " yglob " << global(1)  << " fityglob " << fitpoint(1)
			    << " dzloc " << residual(1)
			    << " X0 " << fitpars[1] << " Y0 " << fitpars[2] 
			    << " derivx R " << lcl_derivativeX[0] << " label " << 1 
			    << " derivx X0 " << lcl_derivativeX[1] << " label " << 2
 			    << " derivx Y0 " << lcl_derivativeX[2] << " label " << 3
			    << " derivx Zslope " << lcl_derivativeX[3] << " label " << 4
			    << " derivx Z0 " << lcl_derivativeX[4] << " label " << 5
			    << " glblderivX alpha " << glbl_derivativeX[0] << " label " << glbl_label[0]
			    << " glblderivX beta " << glbl_derivativeX[1] << " label " << glbl_label[1]
			    << " glblderivX gamma " << glbl_derivativeX[2] << " label " << glbl_label[2]
			    << " glblderivX xtrans " << glbl_derivativeX[3] << " label " << glbl_label[3]
			    << " glblderivX ytrans " << glbl_derivativeX[4] << " label " << glbl_label[4]
			    << " glblderivX ztrans " << glbl_derivativeX[5] << " label " << glbl_label[5]
			    << std::endl;
		}
	    }	  
    
	  if( !isnan(residual(0)) && clus_sigma(0) < 1.0)  // discards crazy clusters
	    { 
	      std::cout << "ckey " << cluskey << " and layer " << layer << " and hitsetkey " << (uint32_t)TrkrDefs::getHitSetKeyFromClusKey(cluskey) << " buffers:" << std::endl;
          AlignmentDefs::printBuffers(0, residual, clus_sigma, lcl_derivativeX, glbl_derivativeX, glbl_label);

	      _mille->mille(AlignmentDefs::NLC, lcl_derivativeX, AlignmentDefs::NGL, glbl_derivativeX, glbl_label, residual(0), errinf*clus_sigma(0));
	    }
	  
	  // provides output that can be grep'ed to make plots of input to mille
	  if(Verbosity() > 1)
	    {
	      //if(layer < 7)
		{
		  std::cout << "Local residualsY: layer " << layer << " phi " << phi * 180 / M_PI << " beta " << beta * 180.90 / M_PI
			    << " dzloc " << residual(1) << " error " << clus_sigma(1) << " inflation_factor " << errinf
			    << " zloc " << zloc << " fitzloc " << fitpoint_local(1)
			    << " zglob " << global(2) << " fitzglob " << fitpoint(2) 
			    << " xglob " << global(0) << " fitxglob " << fitpoint(0)
			    << " yglob " << global(1) << " fityglob " << fitpoint(1)
			    << " dxloc " << residual(0)
			    << " zslope " << fitpars[3] << " Z0 " << fitpars[4]
			    << " derivy R " << lcl_derivativeY[0] << " label " << 1 
			    << " derivy X0 " << lcl_derivativeY[1] << " label " << 2
 			    << " derivy Y0 " << lcl_derivativeY[2] << " label " << 3
			    << " derivy Zslope " << lcl_derivativeY[3] << " label " << 4
			    << " derivy Z0 " << lcl_derivativeY[4] << " label " << 5
			    << " glblderivY alpha " << glbl_derivativeY[0] << " label " << glbl_label[0]
			    << " glblderivY beta " << glbl_derivativeY[1] << " label " << glbl_label[1]
			    << " glblderivY gamma " << glbl_derivativeY[2] << " label " << glbl_label[2]
			    << " glblderivY xtrans " << glbl_derivativeY[3] << " label " << glbl_label[3]
			    << " glblderivY ytrans " << glbl_derivativeY[4] << " label " << glbl_label[4]
			    << " glblderivY ztrans " << glbl_derivativeY[5] << " label " << glbl_label[5]
			    << std::endl;
		}
	    }

	  if(!isnan(residual(1)) && clus_sigma(1) < 1.0 && trkrid != TrkrDefs::inttId)
	    {
	      std::cout << "ckey " << cluskey << " and layer " << layer << " buffers:" << std::endl;
          AlignmentDefs::printBuffers(1, residual, clus_sigma, lcl_derivativeY, glbl_derivativeY, glbl_label);
	      _mille->mille(AlignmentDefs::NLC, lcl_derivativeY, AlignmentDefs::NGL, glbl_derivativeY, glbl_label, residual(1), errinf*clus_sigma(1));
	    }
	}
      tottracks++;
      // close out this track
      _mille->end();
      
    }  // end loop over tracks
  
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 HelicalFitter::get_helix_surface_intersection(Surface surf, std::vector<float>& fitpars, Acts::Vector3 global)
{
  // we want the point where the helix intersects the plane of the surface

  // get the plane of the surface
  Acts::Vector3 sensorCenter      = surf->center(_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
  Acts::Vector3 sensorNormal    = -surf->normal(_tGeometry->geometry().getGeoContext());
  sensorNormal /= sensorNormal.norm();

  // there are analytic solutions for a line-plane intersection.
  // to use this, need to get the vector tangent to the helix near the measurement and a point on it.
  std::pair<Acts::Vector3, Acts::Vector3> line =  get_helix_tangent(fitpars, global);
  Acts::Vector3 pca = line.first;
  Acts::Vector3 tangent = line.second;

  //  std::cout << "   pca: "  << pca(0) << "  " << pca(1) << "  " << pca(2) << "  " << std::endl;

  Acts::Vector3 intersection = get_line_plane_intersection(pca, tangent, sensorCenter, sensorNormal);

  return intersection;
}

Acts::Vector3 HelicalFitter::get_line_plane_intersection(Acts::Vector3 PCA, Acts::Vector3 tangent, Acts::Vector3 sensor_center, Acts::Vector3 sensor_normal)
{
  // get the intersection of the line made by PCA and tangent with the plane of the sensor

  // For a point on the line
  // p = PCA + d * tangent;
  // for a point on the plane
  // (p - sensor_center).sensor_normal = 0

 // The solution is:
  float d = (sensor_center - PCA).dot(sensor_normal) / tangent.dot(sensor_normal);
  Acts::Vector3 intersection = PCA + d * tangent;
  /*
  std::cout << "   intersection: " << intersection(0) << "  "  << intersection(1) << "  "  << intersection(2) << "  " << std::endl;
  std::cout << "        sensor_center: " << sensor_center(0) << "  " << sensor_center(1) << "  " << sensor_center(2) << "  " << std::endl;
  std::cout << "        sensor_normal: " << sensor_normal(0) << "  " << sensor_normal(1) << "  " << sensor_normal(2) << "  " << std::endl;
  */
  return intersection;
}


std::pair<Acts::Vector3, Acts::Vector3> HelicalFitter::get_helix_tangent(const std::vector<float>& fitpars, Acts::Vector3 global)
{
  // no analytic solution for the coordinates of the closest approach of a helix to a point
  // Instead, we get the PCA in x and y to the circle, and the PCA in z to the z vs R line at the R of the PCA 

  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];  
  float zslope = fitpars[3];
  float z0 = fitpars[4];

  Acts::Vector2 pca_circle = get_circle_point_pca(radius, x0, y0, global);

  // The radius of the PCA determines the z position:
  float pca_circle_radius = pca_circle.norm();  // radius of the PCA of the circle to the point
  float pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 pca(pca_circle(0), pca_circle(1), pca_z);

  // now we want a second point on the helix so we can get a local straight line approximation to the track
  // Get the angle of the PCA relative to the fitted circle center
  float angle_pca = atan2(pca_circle(1) - y0, pca_circle(0) - x0);
  // calculate coords of a point at a slightly larger angle
  float d_angle = 0.005;
  float newx = radius * cos(angle_pca + d_angle) + x0;
  float newy = radius * sin(angle_pca + d_angle) + y0;
  float newz = sqrt(newx*newx+newy*newy) * zslope + z0;
  Acts::Vector3 second_point_pca(newx, newy, newz);

  // pca and second_point_pca define a straight line approximation to the track
  Acts::Vector3 tangent = (second_point_pca - pca) /  (second_point_pca - pca).norm();

 // get the PCA of the cluster to that line
  Acts::Vector3 final_pca = getPCALinePoint(global, tangent, pca);

  if(Verbosity() > 2)
    {
      // different method for checking:
      // project the circle PCA vector an additional small amount and find the helix PCA to that point 
      float projection = 0.25;  // cm
      Acts::Vector3 second_point = pca + projection * pca/pca.norm();
      Acts::Vector2 second_point_pca_circle = get_circle_point_pca(radius, x0, y0, second_point);
      float second_point_pca_z = second_point_pca_circle.norm() * zslope + z0;
      Acts::Vector3 second_point_pca2(second_point_pca_circle(0), second_point_pca_circle(1), second_point_pca_z);
      Acts::Vector3 tangent2 = (second_point_pca2 - pca) /  (second_point_pca2 - pca).norm();
      Acts::Vector3 final_pca2 = getPCALinePoint(global, tangent2, pca);
    
      std::cout << " getting tangent at angle_pca: " << angle_pca * 180.0 / M_PI << std::endl 
		<< " pca                      " << pca(0) << "  " << pca(1) << "  " << pca(2) << std::endl
		<< " second_point  " << second_point_pca(0) << "  " << second_point_pca(1) << "  " << second_point_pca(2) << std::endl
		<< " tangent " << tangent(0) << "  " << tangent(1) << "  " << tangent(2) << std::endl	
		<< " final_pca " << final_pca(0) << "  " << final_pca(1) << "  " << final_pca(2) << std::endl	
		<< " second_point2 " << second_point_pca2(0) << "  " << second_point_pca2(1) << "  " << second_point_pca2(2) << std::endl
		<< " tangent2 " << tangent2(0) << "  " << tangent2(1) << "  " << tangent2(2) << std::endl	
		<< " final_pca2 " << final_pca2(0) << "  " << final_pca2(1) << "  " << final_pca2(2) 
		<< std::endl;
    }


  std::pair<Acts::Vector3, Acts::Vector3> line = std::make_pair(final_pca, tangent);

  return line;
}
  
int HelicalFitter::End(PHCompositeNode* )
{
  std::cout << "Tot tracks is " << tottracks<<std::endl;
  _file->cd();
  _residuallayerx->Write();
  _residuallayerz->Write();
  _residuallayerxpull->Write();
  _residuallayerzpull->Write();
  _file->Close();

  // closes output file in destructor
  delete _mille;

  return Fun4AllReturnCodes::EVENT_OK;
}

int  HelicalFitter::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_tpc = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_track_map_tpc)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
    {
      std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
} 

Acts::Vector3 HelicalFitter::get_helix_pca(std::vector<float>& fitpars, Acts::Vector3 global)
{
  return TrackFitUtils::get_helix_pca(fitpars, global);
}


Acts::Vector3 HelicalFitter::getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref)
{
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)   

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 pca = posref + ( (global - posref).dot(tangent) ) * tangent;

  return pca;
}

Acts::Vector2 HelicalFitter::get_circle_point_pca(float radius, float x0, float y0, Acts::Vector3 global)
{
  // get the PCA of a cluster (x,y) position to a circle
  // draw a line from the origin of the circle to the point
  // the intersection of the line with the circle is at the distance radius from the origin along that line 

  Acts::Vector2 origin(x0, y0);
  Acts::Vector2 point(global(0), global(1));

  Acts::Vector2 pca = origin + radius * (point - origin) / (point - origin).norm();

  return pca;
}

float HelicalFitter::convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = _tGeometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89; // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;   // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if(side == 0) zloc = -zloc;
  float z = zloc;  // in cm
 
  return z; 
}

void HelicalFitter::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global)
{
  // make all corrections to global position of TPC cluster
  unsigned int side = TpcDefs::getSide(cluster_key);
  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
  global[2] = z;
  
  // apply distortion corrections
  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
}




// this method to be replaced by calls to TrackFitUtils
void HelicalFitter::getTrackletClusters(TrackSeed *tracklet, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  getTrackletClusterList(tracklet, cluskey_vec);
  // store cluster global positions in a vector
  TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);   
}

void HelicalFitter::getTrackletClusterList(TrackSeed *tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = _cluster_map->findCluster(key);
      if(!cluster)
	{
	  std::cout << "Failed to get cluster with key " << key << std::endl;
	  continue;
	}	  
      
      /// Make a safety check for clusters that couldn't be attached to a surface
      auto surf = _tGeometry->maps().getSurface(key, cluster);
      if(!surf)
	{ continue; }
      
      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}

std::vector<float> HelicalFitter::fitClusters(std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  return TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit
}

Acts::Vector2 HelicalFitter::getClusterError(TrkrCluster *cluster, TrkrDefs::cluskey cluskey, Acts::Vector3& global)
{
  Acts::Vector2 clus_sigma(0,0);

  if(_cluster_version==3)
    {
      clus_sigma(1) = cluster->getZError();
      clus_sigma(0) = cluster->getRPhiError();
    }
  else if(_cluster_version==4)
    {
      double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
      auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster,clusRadius,cluskey);
      float exy2 = para_errors.first * Acts::UnitConstants::cm2;
      float ez2 = para_errors.second * Acts::UnitConstants::cm2;
      clus_sigma(1) = sqrt(ez2);
      clus_sigma(0) = sqrt(exy2);
    }
  else if(_cluster_version == 5)
    {
      double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
      TrkrClusterv5* clusterv5 = dynamic_cast<TrkrClusterv5*>(cluster);
      auto para_errors = _ClusErrPara.get_clusterv5_modified_error(clusterv5,clusRadius,cluskey);
      double phierror = sqrt(para_errors.first);
      double zerror = sqrt(para_errors.second);
      clus_sigma(1) = zerror;
      clus_sigma(0) = phierror;
    }

  return clus_sigma; 
}

// new one
void HelicalFitter::getLocalDerivativesXY(Surface surf, Acts::Vector3 global, const std::vector<float>& fitpars, float lcl_derivativeX[5], float lcl_derivativeY[5], unsigned int layer)
{
  // Calculate the derivatives of the residual wrt the track parameters numerically
  std::vector<float> temp_fitpars;

  std::vector<float> fitpars_delta;
  fitpars_delta.push_back(0.1);   // radius, cm
  fitpars_delta.push_back(0.1);   // X0, cm
  fitpars_delta.push_back(0.1);   // Y0, cm
  fitpars_delta.push_back(0.1);   // zslope, cm
  fitpars_delta.push_back(0.1);   // Z0, cm

  for(unsigned int ip = 0; ip < fitpars.size(); ++ip)
    {
      temp_fitpars.push_back(fitpars[ip]);
    }

  // calculate projX and projY vectors once for the optimum fit parameters
  //  std::pair<Acts::Vector3, Acts::Vector3> tangent = get_helix_tangent(fitpars, global)
  std::pair<Acts::Vector3, Acts::Vector3> tangent = get_helix_tangent(fitpars, global);  // should this be global, not fitpoint?

  Acts::Vector3 projX(0,0,0), projY(0,0,0);
  get_projectionXY(surf, tangent, projX, projY);


  Acts::Vector3 intersection = get_helix_surface_intersection(surf, temp_fitpars, global);

  // loop over the track fit parameters
  for(unsigned int ip = 0; ip < fitpars.size(); ++ip)
    {
      temp_fitpars[ip] += fitpars_delta[ip];

      Acts::Vector3 temp_intersection = get_helix_surface_intersection(surf, temp_fitpars, global);
      Acts::Vector3 intersection_delta = temp_intersection - intersection;
      if(Verbosity() > 1)
	{
	  std::cout << "Layer " << layer << " local parameter " << ip << ":" << std::endl; 
	  std::cout << " intersection " << intersection(0) << "  " << intersection(1) << "  " << intersection(2) << std::endl;
	  std::cout << " temp_intersection " << temp_intersection(0) << "  "<< temp_intersection(1) << "  "<< temp_intersection(2)<< std::endl;
	  std::cout << " intersection_delta " << intersection_delta(0) << "  " << intersection_delta(1) << "  " << intersection_delta(2) << std::endl;
	}


      // convert to delta-intersection / delta-parameter
      intersection_delta /= fitpars_delta[ip];

      if(Verbosity() > 1)
	{std::cout << " intersection_delta / delta_p " << intersection_delta(0) << "  " << intersection_delta(1) << "  " << intersection_delta(2) << std::endl;}

      // calculate the change in residual for X and Y 
      lcl_derivativeX[ip] = intersection_delta.dot(projX);
      lcl_derivativeY[ip] = intersection_delta.dot(projY);
      if(Verbosity() > 1)
	{std::cout << " ip " << ip << "  derivativeX " << lcl_derivativeX[ip] << "  " << " derivativeY " << lcl_derivativeY[ip] << std::endl;}

      temp_fitpars[ip] = fitpars[ip];
    }
}

void HelicalFitter::getGlobalDerivativesXY(Surface surf, Acts::Vector3 global, Acts::Vector3 fitpoint, const std::vector<float>& fitpars, float glbl_derivativeX[6], float glbl_derivativeY[6], unsigned int layer)
{
  Acts::Vector3 unitx(1, 0, 0);
  Acts::Vector3 unity(0, 1, 0);
  Acts::Vector3 unitz(0, 0, 1);

  // calculate projX and projY vectors once for the optimum fit parameters
  std::pair<Acts::Vector3, Acts::Vector3> tangent = get_helix_tangent(fitpars, global);  // should this be global, not fitpoint?

  Acts::Vector3 projX(0,0,0), projY(0,0,0);
  get_projectionXY(surf, tangent, projX, projY);

  // translations

  glbl_derivativeX[3] = unitx.dot(projX);
  glbl_derivativeX[4] = unity.dot(projX);
  glbl_derivativeX[5] = unitz.dot(projX);

  glbl_derivativeY[3] = unitx.dot(projY);
  glbl_derivativeY[4] = unity.dot(projY);
  glbl_derivativeY[5] = unitz.dot(projY);

  // rotations

  // need center of sensor to intersection point
  Acts::Vector3 sensorCenter      = surf->center(_tGeometry->geometry().getGeoContext()) / Acts::UnitConstants::cm;  // convert to cm
  Acts::Vector3 OM = fitpoint - sensorCenter;

  glbl_derivativeX[0] = (unitx.cross(OM)).dot(projX);
  glbl_derivativeX[1] = (unity.cross(OM)).dot(projX);
  glbl_derivativeX[2] = (unitz.cross(OM)).dot(projX);

  glbl_derivativeY[0] = (unitx.cross(OM)).dot(projY);
  glbl_derivativeY[1] = (unity.cross(OM)).dot(projY);
  glbl_derivativeY[2] = (unitz.cross(OM)).dot(projY);

  if(Verbosity() > 1)
    {
      std::cout << " glbl_derivativesX for layer " << layer << std::endl;
      for(unsigned int i = 0; i < 6; ++i)
	{
	  std::cout << " i " << i << " glbl_derivative " << glbl_derivativeX[i] << std::endl;
	}
      
      std::cout << " glbl_derivativesY for layer " << layer << std::endl;
      for(unsigned int i = 0; i < 6; ++i)
	{
	  std::cout << " i " << i << " glbl_derivative " << glbl_derivativeY[i] << std::endl;
	}
    }


}
void HelicalFitter::get_projectionXY(Surface surf, std::pair<Acts::Vector3, Acts::Vector3> tangent, Acts::Vector3& projX, Acts::Vector3& projY)
{
  // we only need the direction part of the tangent
  Acts::Vector3 tanvec = tangent.second;

  // get the plane of the surface
  Acts::Vector3 sensorCenter      = surf->center(_tGeometry->geometry().getGeoContext()) / Acts::UnitConstants::cm;  // convert to cm
  // sensorNormal is the Z vector
  Acts::Vector3 Z = -surf->normal(_tGeometry->geometry().getGeoContext()) / Acts::UnitConstants::cm; 

  // get surface X and Y unit vectors in global frame
  // transform Xlocal = 1.0 to global, subtract the surface center, normalize to 1
  Acts::Vector3 xloc(1.0,0.0,0.0);
  Acts::Vector3 xglob =  surf->transform(_tGeometry->geometry().getGeoContext()) * (xloc *  Acts::UnitConstants::cm);
  xglob /=  Acts::UnitConstants::cm;
  Acts::Vector3 yloc(0.0,1.0,0.0);
  Acts::Vector3 yglob =  surf->transform(_tGeometry->geometry().getGeoContext()) * (yloc *  Acts::UnitConstants::cm);
  yglob /=  Acts::UnitConstants::cm;

  Acts::Vector3 X = (xglob-sensorCenter) / (xglob-sensorCenter).norm();
  Acts::Vector3 Y = (yglob-sensorCenter) / (yglob-sensorCenter).norm();

  projX = X - (tanvec.dot(X) / tanvec.dot(Z)) * Z;
  projY = Y - (tanvec.dot(Y) / tanvec.dot(Z)) * Z;

  return;
}

unsigned int HelicalFitter::addSiliconClusters(std::vector<float>& fitpars, std::vector<Acts::Vector3>& global_vec,  std::vector<TrkrDefs::cluskey>& cluskey_vec)
{

  return TrackFitUtils::addSiliconClusters(fitpars, dca_cut, _tGeometry, _cluster_map, global_vec, cluskey_vec);
}

 bool HelicalFitter::is_layer_fixed(unsigned int layer)
 {
  bool ret = false;
  auto it = fixed_layers.find(layer);
  if(it != fixed_layers.end()) 
    ret = true;

  return ret;
 }

 void HelicalFitter::set_layer_fixed(unsigned int layer)
 {
   fixed_layers.insert(layer);
 }

 
bool HelicalFitter::is_layer_param_fixed(unsigned int layer, unsigned int param)
 {
  bool ret = false;
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  auto it = fixed_layer_params.find(pair);
  if(it != fixed_layer_params.end()) 
    ret = true;

  return ret;
 }

void HelicalFitter::set_layer_param_fixed(unsigned int layer, unsigned int param)
 {
   std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
   fixed_layer_params.insert(pair);
 }

void HelicalFitter::set_tpc_sector_fixed(unsigned int region, unsigned int sector, unsigned int side)
 {
   // make a combined subsector index
   unsigned int subsector = region * 24 + side * 12 + sector;
   fixed_sectors.insert(subsector);
 }

bool HelicalFitter::is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side)
 {
   bool ret = false;
   unsigned int region = AlignmentDefs::getTpcRegion(layer);
   unsigned int subsector = region * 24 + side * 12 + sector;
   auto it = fixed_sectors.find(subsector);
   if(it != fixed_sectors.end()) 
     ret = true;
  
   return ret;
 }

void HelicalFitter::correctTpcGlobalPositions(std::vector<Acts::Vector3> global_vec,  std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  for(unsigned int iclus=0;iclus<cluskey_vec.size();++iclus)
    {
      auto cluskey = cluskey_vec[iclus];
      auto global = global_vec[iclus];
      const unsigned int trkrId = TrkrDefs::getTrkrId(cluskey);	  
      if(trkrId == TrkrDefs::tpcId) 
	{  
	  // have to add corrections for TPC clusters after transformation to global
	  int crossing = 0;  // for now
	  makeTpcGlobalCorrections(cluskey, crossing, global); 
	  global_vec[iclus] = global;
	}
    }
}

