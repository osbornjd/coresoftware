#include "PHCircleFit.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>

PHCircleFit::PHCircleFit(const std::string& name)
  : SubsysReco(name)
{
}


int PHCircleFit::InitRun(PHCompositeNode *topNode)
{
  int returnval = getNodes(topNode);

  return returnval;
}

int PHCircleFit::process_event(PHCompositeNode *topNode)
{

  for(const auto& [key, track] : *m_trackMap)
    {
      if(Verbosity() > 3)
	{ track->identify(); }
      
      Acts::Vector3D position, momentum;
      int charge;

      circleFitTrack(track, position, momentum, charge);

      track->set_x(position(0));
      track->set_y(position(1));
      track->set_z(position(2));
      track->set_px(momentum(0));
      track->set_py(momentum(1));
      track->set_pz(momentum(2));
      track->set_charge(charge);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHCircleFit::circleFitTrack(SvtxTrack *track, Acts::Vector3D& position,
				 Acts::Vector3D& momentum, int& charge)
{

  auto clusters = getClusters(track);

  double R, X0, Y0;
  circleFitByTaubin(clusters, R, X0, Y0);
  
  double x, y;
  findRoot(R, X0, Y0, x, y,
	   track->get_vertex_id());

  position(0) = x;
  position(1) = y;
  
  charge = getCharge(clusters, atan2(Y0,X0));
  
  /// Now determine the line tangent to the circle at this point to get phi
  /// The slope of the line connecting the circle center and PCA is 
  /// m = (y0-y)/(x0-x). So the perpendicular slope (i.e. phi) is then -1/m
  /// For some reason the phi value comes back from atan2 off by 
  /// a factor of pi for positive charged tracks, hence the check
  
  double phi = atan2(-1 * (X0-x), Y0-y);
  if(charge > 0)
    {
      phi += M_PI;
      if(phi > M_PI) 
	phi -= 2. * M_PI;
    }
 
  if(Verbosity() > 2)
    std::cout << "track seed phi : " << phi <<  std::endl;

  double m, B;
  
  /// m is slope as a function of radius, B is z intercept (vertex)
  
  lineFit(clusters, m, B);

  /// Need to evaluate position at r intercept given by x,y positions
  position(2) = evaluateZVertex(track->get_vertex_id(), m, B, x, y);
  
  double theta = atan(1./m);

  /// normalize to 0 < theta < pi
  if(theta < 0)
    theta += M_PI;

  if(Verbosity() > 2)
    { std::cout << "Track seed theta: " << theta << std::endl; }
 
  /// 0.3 conversion factor, 1.4=B field, 100 convert R from cm to m
  /// Get a very rough estimate of p
  float pt = 0.3 * 1.4 * R / 100.;
  float eta = -log(tan(theta/2.));
  float p = pt * cosh(eta);

  /// The only thing that is really needed for the propagation
  /// is the direction
  momentum(0) = p * sin(theta) * cos(phi);
  momentum(1) = p * sin(theta) * sin(phi);
  momentum(2) = p * cos(theta);
  
  if(Verbosity() > 2)
    {  
      std::cout << "Phi,eta : " << phi << " , " << eta << std::endl;
      std::cout << "Momentum vector estimate: " << momentum.transpose() 
		<< std::endl; 
      std::cout << "Position estimate: " << position.transpose() 
		<< std::endl;
    }
 

}

double PHCircleFit::evaluateZVertex(const int vtxid, double& m, 
				    double& B, double& x, double& y)
{
  auto vertex = m_vertexMap->get(vtxid);
  double vz = vertex->get_z();
  double r = sqrt(x*x + y*y);

  /// Because R can only be positive, we can have two z solutions
  double z1 = m * r + B;
  double z2 = -z1;

  /// determine which is closer to the z vertex position
  double dz1 = fabs(vz - z1);
  double dz2 = fabs(vz - z2);
  
  std::cout << "vz " << vz << "z1 z2 " << z1 << ", " << z2 << " dz1,dz2" 
	    << dz1 << ", " << dz2 << std::endl;
  if(dz1 < dz2)
    return z1;

  return z2;

}

void PHCircleFit::lineFit(const std::vector<TrkrCluster*>& clusters,
			  double& A, double& B)
{

 // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
  double xsum = 0,x2sum = 0,ysum = 0,xysum = 0;
  int clussize = 0;
  for(auto& cluster : clusters)
    {
      /// Intt has bad z resolution. Screws up fit
      if(TrkrDefs::getTrkrId(cluster->getClusKey()) == TrkrDefs::inttId)
	{ continue; }
      
      clussize++;
      double z = cluster->getZ();
      double r = sqrt(pow(cluster->getX(),2) + pow(cluster->getY(), 2));
      
      xsum=xsum+r;               // calculate sigma(xi)
      ysum=ysum+z;               // calculate sigma(yi)
      x2sum=x2sum+pow(r,2);      // calculate sigma(x^2i)
      xysum=xysum+r*z;           // calculate sigma(xi*yi)
    }
  
  /// calculate slope
  A = (clussize*xysum-xsum*ysum) / (clussize*x2sum-xsum*xsum);

  /// calculate intercept
  B = (x2sum*ysum-xsum*xysum) / (x2sum*clussize-xsum*xsum);
  
  if(Verbosity() > 2)
    {
      for(auto& cluster : clusters)
	{
	  
	  double r = sqrt(pow(cluster->getX(),2) + pow(cluster->getY(),2));
	  double zfit = A * r + B;
	  std::cout << "r " << r << ", clus z " << cluster->getZ() 
		    << ", z fit " << zfit << std::endl;
	}

    }
  

}

int PHCircleFit::getCharge(const std::vector<TrkrCluster*>& clusters,
			   const double& circPhi)
{
  
  /**
   * If the circle center phi is positioned clockwise to the seed phi, 
   * the seed is positively charged. If the circle center phi is positioned
   * counter clockwise, the seed is negatively charged
   */

  int charge = 0;
  
  /// Get a crude estimate of the seed phi by taking the average of the
  /// measurements
  double trackPhi = 0;
  for(auto& clus : clusters)
    {
      double clusPhi = atan2(clus->getY(), clus->getX());

      /// if it is close to the periodic boundary normalize to 
      /// two pi to avoid -pi and pi issues
      if(fabs(fabs(clusPhi) - M_PI) < 0.2)
	clusPhi = normPhi2Pi(clusPhi);
      trackPhi += clusPhi;
    }

  trackPhi /= clusters.size();

  /// normalize back
  if(trackPhi > M_PI)
    trackPhi -= 2. * M_PI;

  float quadrants[5] = {-M_PI,-M_PI / 2., 0, M_PI/2., M_PI};
  int quadrant = -1;
  for(int i=0; i<4; i++)
    {
      if(trackPhi > quadrants[i] && trackPhi <= quadrants[i+1])
	{
	  quadrant = i;
	  break;
	}
    }

  if(quadrant == -1)
    std::cout << "quadrant was not set... shouldn't be possible"
	      << std::endl;

  if(quadrant == 1 or quadrant == 2)
    {
      if(circPhi > trackPhi)
	charge = -1;
      else
	charge = 1;
    }
  else
    {
      /// Shift the periodic boundary to make quadrants 0 and 3 away
      /// from boundary
      double normTrackPhi = normPhi2Pi(trackPhi);
      double normCircPhi = normPhi2Pi(circPhi);
  
      if(normCircPhi > normTrackPhi)
	charge = -1;
      else
	charge = 1;
    }

  if(Verbosity() > 1)
    std::cout << "Track seed charge determined to be " 
	      << charge << " in quadrant " << quadrant << std::endl;

  return charge;

}

void PHCircleFit::findRoot(double& R, double& X0, double& Y0,
			   double& x, double& y, const int& vertexId)
{
    /**
   * We need to determine the closest point on the circle to the vertex
   * since we can't assume that the track originates from the vertex
   * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to 
   * minimize d = sqrt((vx-x)^2+(vy-y)^2), the distance between the 
   * vertex and some (currently, unknown) point on the circle x,y.
   * 
   * Solving the circle eqn for x and substituting into d gives an eqn for
   * y. Taking the derivative and setting equal to 0 gives the following 
   * two solutions for an arbitrary vertex (vx,vy). This condenses to 
   * the solution in PHActsSiliconSeeding for vx=vy=0.
   */
  
  auto vertex = m_vertexMap->get(vertexId);
  
  float vx = 0, vy = 0;

  if(vertex != nullptr)
    {
      vx = vertex->get_x();
      vy = vertex->get_y();
      if(Verbosity() > 2)
	{ std::cout << "Vertex values are " << vx << ", " << vy << std::endl; }
    }

  double a = 2*X0*vx-vx*vx-X0*X0-vy*vy+2*vy*Y0-Y0*Y0;
  double b = 2*Y0*Y0*Y0 - 4*Y0*Y0*vy + 2*Y0*vy*vy+2*X0*X0*Y0-4*X0*Y0*vx+2*Y0*vx*vx;
  double c = 2*Y0*Y0*Y0*vy-Y0*Y0*Y0*Y0-Y0*Y0*vy*vy+R*R*Y0*Y0-2*R*R*Y0*vy+R*R*vy*vy-X0*X0*Y0*Y0+2*X0*Y0*Y0*vx-Y0*Y0*vx*vx;

  double y1 = (-1 * b - sqrt(b*b - 4*a*c)) / (2.*a);
  double y2 = (-1 * b + sqrt(b*b - 4*a*c)) / (2.*a);
  double x1 = sqrt(pow(R, 2) - pow(y1 - Y0, 2)) + X0;
  double x2 = -sqrt(pow(R, 2) - pow(y2 - Y0, 2)) + X0;

  if(fabs(x1) < fabs(x2))
    x = x1;
  else
    x = x2;

  if(fabs(y1) < fabs(y2))
    y = y1;
  else
    y = y2;

  if(Verbosity() > 1)
    std::cout << "x1 and x2 : " << x1 << ", " << x2 << std::endl
	      << "y1 and y2 : " << y1 << ", " << y2 << std::endl;

  if(Verbosity() > 1)
    {
      std::cout << "Minimum x and y positions " << x << ",  " 
		<< y << std::endl;
    }


}

void PHCircleFit::circleFitByTaubin(std::vector<TrkrCluster*>& clusters,
				    double& R, double& X0, double& Y0)
{
   /**  
   *   Circle fit to a given set of data points (in 2D)
   *   This is an algebraic fit, due to Taubin, based on the journal article
   *   G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
   *               Space Curves Defined By Implicit Equations, With 
   *               Applications To Edge And Range Image Segmentation",
   *               IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
   *  It works well whether data points are sampled along an entire circle 
   *  or along a small arc. 
   *  It still has a small bias and its statistical accuracy is slightly lower 
   *  than that of the geometric fit (minimizing geometric distances),
   *  It provides a very good initial guess for a subsequent geometric fit. 
   *    Nikolai Chernov  (September 2012)
   */
  
  int iter, IterMAX=99;
  
  double Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
  double A0, A1, A2, A22, A3, A33;
  double x, y;
  double DET, Xcenter, Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;
  
  for(auto& clus : clusters)
    {
      meanX += clus->getX();
      meanY += clus->getY();
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

  for(auto& clus : clusters)
    {
      double Xi = clus->getX() - meanX;
      double Yi = clus->getY() - meanY;
      double Zi = Xi * Xi + Yi * Yi;

      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
    }

  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;

  Mz = Mxx + Myy;
  Cov_xy = Mxx * Myy - Mxy * Mxy;
  Var_z = Mzz - Mz * Mz;
  A3 = 4 * Mz;
  A2 = -3 * Mz * Mz - Mzz;
  A1 = Var_z * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
  A0 = Mxz * (Mxz * Myy - Myz * Mxy) + 
    Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;

  for (x=0., y=A0, iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x * (A22 + A33 * x);
      double xnew = x - y / Dy;
      if ((xnew == x)||(!std::isfinite(xnew))) break;
      double ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  assembling the output
  
  X0 = Xcenter + meanX;
  Y0 = Ycenter + meanY;
  R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);


}

std::vector<TrkrCluster*> PHCircleFit::getClusters(SvtxTrack *track)
{
  std::vector<TrkrCluster*> clusters;

  for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
    {
      auto cluster = m_clusterContainer->findCluster(*iter);
      /// Skip whatever layers are requested to skip
      auto layer = TrkrDefs::getLayer(cluster->getClusKey());
  
      if(layer < m_minLayer or layer > m_maxLayer)
	{ continue; }
      
      clusters.push_back(cluster);
    }

  return clusters;
}



int PHCircleFit::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


double PHCircleFit::normPhi2Pi(const double& phi)
{
  double returnPhi = phi;
  if(returnPhi < 0)
    returnPhi += 2 * M_PI;
  return returnPhi;
}


int PHCircleFit::getNodes(PHCompositeNode *topNode)
{
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "Can't do circle fit without vertex estimate, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "Can't do circle fit without tracks, exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if(!m_clusterContainer)
    {
      std::cout << PHWHERE << "Can't do circle fit without clusters, exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}
