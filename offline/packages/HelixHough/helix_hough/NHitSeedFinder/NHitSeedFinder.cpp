#include "NHitSeedFinder.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <Eigen/LU>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

NHitSeedFinder::NHitSeedFinder(vector<float>& detrad, unsigned int n_phi,
                               unsigned int n_d, unsigned int n_k,
                               unsigned int n_dzdl, unsigned int n_z0,
                               HelixResolution& min_resolution, 
                               HelixResolution& max_resolution, HelixRange& range) : 
  HelixHough(n_phi, n_d, n_k, n_dzdl, n_z0, min_resolution, max_resolution, range),
  using_vertex(false),
  vertex_sigma_xy(0.002),
  vertex_sigma_z(0.005),
  chi2_cut(3.),
  seed_mode(4)
{
  for(unsigned int i=0;i<detrad.size();++i)
  {
    detector_radii.push_back(detrad[i]);
  }
}

void NHitSeedFinder::findTracks(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks)
{
	if(seed_mode == 4)
		find4Tracks(hits,tracks);
	if(seed_mode == 5)
		find5Tracks(hits,tracks);
	if(seed_mode == 6)
		find6Tracks(hits,tracks);
}


void NHitSeedFinder::find4Tracks(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks)
{
  // loop over the provided hits
  for(unsigned int i=0;i<hits.size();i++)
  {
    // if this index is not in our map
    if( phis.find(hits[i].index) == phis.end() )
    {
      // calculate the phi now and insert into the map
      phis.insert(make_pair(hits[i].index,atan2(hits[i].y,hits[i].x)));
    }
  }
	
  vector<double> chi2_hits;
  SimpleTrack3D temp_track;
  temp_track.hits.resize(4, hits[0]);
	
  for(unsigned int i1=0;i1<hits.size();++i1)
    {
      temp_track.hits[0] = hits[i1];
      
      for(unsigned int i2=(i1+1);i2<hits.size();++i2)
				{
					if( hits[i2].layer == hits[i1].layer ) continue;
					if( fabs( phis[hits[i2].index] - phis[hits[i1].index] ) > M_PI_2 ) continue; 
					
					temp_track.hits[1] = hits[i2];
					
					for(unsigned int i3=(i2+1);i3<hits.size();++i3)
						{
							if( hits[i3].layer == hits[i2].layer ) continue;
							if( hits[i3].layer == hits[i1].layer ) continue;
							if( fabs( phis[hits[i3].index] - phis[hits[i2].index] ) > M_PI_2 ) continue; 
							if( fabs( phis[hits[i3].index] - phis[hits[i1].index] ) > M_PI_2 ) continue; 
							
							temp_track.hits[2] = hits[i3];
							
							for(unsigned int i4=(i3+1);i4<hits.size();++i4)
								{
									if( hits[i4].layer == hits[i3].layer ) continue;
									if( hits[i4].layer == hits[i2].layer ) continue;
									if( hits[i4].layer == hits[i1].layer ) continue;
									if( fabs( phis[hits[i4].index] - phis[hits[i3].index] ) > M_PI_2 ) continue; 
									if( fabs( phis[hits[i4].index] - phis[hits[i2].index] ) > M_PI_2 ) continue; 
									if( fabs( phis[hits[i4].index] - phis[hits[i1].index] ) > M_PI_2 ) continue; 
									
									temp_track.hits[3] = hits[i4];
									
									// build track
									vector<unsigned int> tempcomb;
									tempcomb.assign(4,0);
									tempcomb[0] = temp_track.hits[0].index;
									tempcomb[1] = temp_track.hits[1].index;
									tempcomb[2] = temp_track.hits[2].index;
									tempcomb[3] = temp_track.hits[3].index;
									
									// check this combo has already been seen
									sort(tempcomb.begin(), tempcomb.end());
									set<vector<unsigned int> >::iterator it = combos.find(tempcomb);
									if(it != combos.end()) continue;
									
									// add to list of checked combinations
									combos.insert(tempcomb);
									
									// fit track
									double chi2 = fitTrack(temp_track, chi2_hits);
									if(chi2 > chi2_cut) continue;
									
									// add track
									tracks.push_back(temp_track);
								}
						}
				}
    }
}


void NHitSeedFinder::find5Tracks(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks)
{//cout<<"find6Tracks Entered with "<<hits.size()<<" hits"<<endl;                                                     
  // loop over the provided hits                                                                                      
  for(unsigned int i=0;i<hits.size();i++)
    {
      // if this index is not in our map                                                                              
      if( phis.find(hits[i].index) == phis.end() )
        {
          // calculate the phi now and insert into the map                                                            
          phis.insert(make_pair(hits[i].index,atan2(hits[i].y,hits[i].x)));
        }
    }

  vector<double> chi2_hits;
  SimpleTrack3D temp_track;
  temp_track.hits.resize(5, hits[0]);

  for(unsigned int i1=0;i1<hits.size();++i1)
    {
      temp_track.hits[0] = hits[i1];

      for(unsigned int i2=(i1+1);i2<hits.size();++i2)
        {
          if( hits[i2].layer == hits[i1].layer ) continue;
					/*              if( fabs( phis[hits[i3].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;            
													if( fabs( phis[hits[i3].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

					temp_track.hits[1] = hits[i2];

					for(unsigned int i3=(i2+1);i3<hits.size();++i3)
						{
							if( hits[i3].layer == hits[i2].layer ) continue;
							if( hits[i3].layer == hits[i1].layer ) continue;
							/*                  if( fabs( phis[hits[i4].index] - phis[hits[i3].index] ) > M_PI_2 ) continue;    
                                  if( fabs( phis[hits[i4].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;                        
                                  if( fabs( phis[hits[i4].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

							temp_track.hits[2] = hits[i3];

							for(unsigned int i4=(i3+1);i4<hits.size();++i4)
								{
									if( hits[i4].layer == hits[i3].layer ) continue;
									if( hits[i4].layer == hits[i2].layer ) continue;
									if( hits[i4].layer == hits[i1].layer ) continue;
									/*                      if( fabs( phis[hits[i5].index] - phis[hits[i4].index] ) > M_PI_2 ) continue;                                                                                                                 
                                          if( fabs( phis[hits[i5].index] - phis[hits[i3].index] ) > M_PI_2 ) continue;                    
                                          if( fabs( phis[hits[i5].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;                    
                                          if( fabs( phis[hits[i5].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

									temp_track.hits[3] = hits[i4];
									for(unsigned int i5=(i4+1);i5<hits.size();++i5)
										{
											if( hits[i5].layer == hits[i4].layer ) continue;
											if( hits[i5].layer == hits[i3].layer ) continue;
											if( hits[i5].layer == hits[i2].layer ) continue;
											if( hits[i5].layer == hits[i1].layer ) continue;
											/*                      if( fabs( phis[hits[i5].index] - phis[hits[i4].index] ) > M_PI_2 ) continue;                                                           if( fabs( phis[hits[i5].index] - phis[hits[i3].index] ) > M_PI_2 ) continue;                                                                       if( fabs( phis[hits[i5].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;                                                            if( fabs( phis[hits[i5].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/
											temp_track.hits[4]= hits[i5];
											// build track                                                                           
                                                                                                                      
                      vector<unsigned int> tempcomb;
                      tempcomb.assign(5,0);
                      tempcomb[0] = temp_track.hits[0].index;
                      tempcomb[1] = temp_track.hits[1].index;
                      tempcomb[2] = temp_track.hits[2].index;
                      tempcomb[3] = temp_track.hits[3].index;
                      tempcomb[4] = temp_track.hits[4].index;

                      // check this combo has already been seen                                                 
                                                                                                                      
                      sort(tempcomb.begin(), tempcomb.end());
                      set<vector<unsigned int> >::iterator it = combos.find(tempcomb);
                      if(it != combos.end()) continue;
													
                      // add to list of checked combinations                                          
														
                      combos.insert(tempcomb);
														
                      // fit track                                                           
															
                      double chi2 = fitTrack(temp_track, chi2_hits);
                      //                          cout <<"chi2 of seed track: "<<chi2<<endl;                      
                      if(chi2 > chi2_cut) continue;
															
                      // add track                                                                              
																tracks.push_back(temp_track);
										}
								}
						}
				}
		}
}



void NHitSeedFinder::find6Tracks(vector<SimpleHit3D>& hits, vector<SimpleTrack3D>& tracks)
{
//   cout << "find6Tracks Entered with " << hits.size()
//        << " hits, ncombos: " << combos.size() << endl;

  // loop over the provided hits
  for(unsigned int i=0;i<hits.size();i++)
    {
      // if this index is not in our map
      if( phis.find(hits[i].index) == phis.end() )
				{
					// calculate the phi now and insert into the map
					phis.insert(make_pair(hits[i].index,atan2(hits[i].y,hits[i].x)));
				}
    }

  vector<double> chi2_hits;
  SimpleTrack3D temp_track;
  temp_track.hits.resize(6, hits[0]);

  for(unsigned int i1=0;i1<hits.size();++i1)
    {
      temp_track.hits[0] = hits[i1];

      for(unsigned int i2=(i1+1);i2<hits.size();++i2)
				{
					if( hits[i2].layer == hits[i1].layer ) continue;

					//if( fabs( phis[hits[i2].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;

					temp_track.hits[1] = hits[i2];

					for(unsigned int i3=(i2+1);i3<hits.size();++i3)
						{
							if( hits[i3].layer == hits[i2].layer ) continue;
							if( hits[i3].layer == hits[i1].layer ) continue;
							/*							if( fabs( phis[hits[i3].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;
								if( fabs( phis[hits[i3].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

							temp_track.hits[2] = hits[i3];

							for(unsigned int i4=(i3+1);i4<hits.size();++i4)
								{
									if( hits[i4].layer == hits[i3].layer ) continue;
									if( hits[i4].layer == hits[i2].layer ) continue;
									if( hits[i4].layer == hits[i1].layer ) continue;
									/*									if( fabs( phis[hits[i4].index] - phis[hits[i3].index] ) > M_PI_2 ) continue;
									if( fabs( phis[hits[i4].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;
									if( fabs( phis[hits[i4].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

									temp_track.hits[3] = hits[i4];

									for(unsigned int i5=(i4+1);i5<hits.size();++i5)
										{
											if( hits[i5].layer == hits[i4].layer ) continue;
											if( hits[i5].layer == hits[i3].layer ) continue;
											if( hits[i5].layer == hits[i2].layer ) continue;
											if( hits[i5].layer == hits[i1].layer ) continue;
											/*                      if( fabs( phis[hits[i5].index] - phis[hits[i4].index] ) > M_PI_2 ) continue;
											if( fabs( phis[hits[i5].index] - phis[hits[i3].index] ) > M_PI_2 ) continue;
											if( fabs( phis[hits[i5].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;
											if( fabs( phis[hits[i5].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

											temp_track.hits[4] = hits[i5];

											for(unsigned int i6=(i5+1);i6<hits.size();++i6)
												{
													if( hits[i6].layer == hits[i5].layer ) continue;
													if( hits[i6].layer == hits[i4].layer ) continue;
													if( hits[i6].layer == hits[i3].layer ) continue;
													if( hits[i6].layer == hits[i2].layer ) continue;
													if( hits[i6].layer == hits[i1].layer ) continue;
													/*													if( fabs( phis[hits[i6].index] - phis[hits[i5].index] ) > M_PI_2 ) continue;
													if( fabs( phis[hits[i6].index] - phis[hits[i4].index] ) > M_PI_2 ) continue;
													if( fabs( phis[hits[i6].index] - phis[hits[i3].index] ) > M_PI_2 ) continue;
													if( fabs( phis[hits[i6].index] - phis[hits[i2].index] ) > M_PI_2 ) continue;
													if( fabs( phis[hits[i6].index] - phis[hits[i1].index] ) > M_PI_2 ) continue;*/

													temp_track.hits[5] = hits[i6];
													
													// build track
													vector<unsigned int> tempcomb;
													tempcomb.assign(6,0);
													tempcomb[0] = temp_track.hits[0].index;
													tempcomb[1] = temp_track.hits[1].index;
													tempcomb[2] = temp_track.hits[2].index;
													tempcomb[3] = temp_track.hits[3].index;
													tempcomb[4] = temp_track.hits[4].index;
                          tempcomb[5] = temp_track.hits[5].index;
													
													// check this combo has already been seen
													sort(tempcomb.begin(), tempcomb.end());
													set<vector<unsigned int> >::iterator it = combos.find(tempcomb);
													if(it != combos.end()) {
                            //cout << "combo already used!" << endl;
                            continue;
													}

													// add to list of checked combinations
													combos.insert(tempcomb);

													// fit track
													double chi2 = fitTrack(temp_track, chi2_hits);
                          //cout <<"chi2 of seed track: "<<chi2<<endl;
													if(chi2 > chi2_cut) continue;

													
													// add track
													tracks.push_back(temp_track);
												}
										}
								}
						}
				}
		}
}


void NHitSeedFinder::finalize(vector<SimpleTrack3D>& input, vector<SimpleTrack3D>& output)
{
  unsigned int nt = input.size();
  for(unsigned int i=0;i<nt;++i)
  {
    output.push_back(input[i]);
  }
}


double NHitSeedFinder::fitTrack(SimpleTrack3D& track, vector<double>& chi2_hit)
{
  if(using_vertex == true)
  {
    track.hits.push_back(SimpleHit3D(0.,0., 0.,0., 0.,0., 0, 0));
  }
  
  chi2_hit.clear();
  chi2_hit.resize(track.hits.size(), 0.);
  
  MatrixXf y = MatrixXf::Zero(track.hits.size(), 1);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    y(i, 0) = ( pow(track.hits[i].x,2) + pow(track.hits[i].y,2) );
    if((using_vertex==true ) && (i == (track.hits.size() - 1))){y(i, 0) /= vertex_sigma_xy;}
    else{y(i, 0) /= layer_xy_resolution[track.hits[i].layer];}
  }
  
  MatrixXf X = MatrixXf::Zero(track.hits.size(), 3);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    X(i, 0) = track.hits[i].x;
    X(i, 1) = track.hits[i].y;
    X(i, 2) = -1.;
    if((using_vertex==true ) && (i == (track.hits.size() - 1)))
    {
      X(i, 0) /= vertex_sigma_xy;
      X(i, 1) /= vertex_sigma_xy;
      X(i, 2) /= vertex_sigma_xy;
    }
    else
    {
      X(i, 0) /= layer_xy_resolution[track.hits[i].layer];
      X(i, 1) /= layer_xy_resolution[track.hits[i].layer];
      X(i, 2) /= layer_xy_resolution[track.hits[i].layer];
    }
  }
  
  MatrixXf Xt = X.transpose();
  
  MatrixXf prod = Xt*X;
  MatrixXf inv = prod.fullPivLu().inverse();
  
  MatrixXf beta = inv*Xt*y;
  
  float cx = beta(0,0)*0.5;
  float cy = beta(1,0)*0.5;
  float r = sqrt(cx*cx + cy*cy - beta(2,0));
  
  float phi = atan2(cy, cx);
  float d = sqrt(cx*cx + cy*cy) - r;
  float k = 1./r;
  
  MatrixXf diff = y - (X*beta);
  MatrixXf chi2 = (diff.transpose())*diff;
  
  float dx = d*cos(phi);
  float dy = d*sin(phi);
  
  MatrixXf y2 = MatrixXf::Zero(track.hits.size(), 1);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    y2(i,0) = track.hits[i].z;
    if((using_vertex==true ) && (i == (track.hits.size() - 1))){y2(i, 0) /= vertex_sigma_z;}
    else{y2(i, 0) /= layer_z_resolution[track.hits[i].layer];}
  }
  
  MatrixXf X2 = MatrixXf::Zero(track.hits.size(), 2);
  for(unsigned int i=0;i<track.hits.size();i++)
  {
    float D = sqrt( pow(dx - track.hits[i].x, 2) + pow(dy - track.hits[i].y,2));
    float s = 0.0;
    
    if(0.5*k*D > 0.1)
    {
      float v = 0.5*k*D;
      if(v >= 0.999999){v = 0.999999;}
      s = 2.*asin(v)/k;
    }
    else
    {
      float temp1 = k*D*0.5;temp1*=temp1;
      float temp2 = D*0.5;
      s += 2.*temp2;
      temp2*=temp1;
      s += temp2/3.;
      temp2*=temp1;
      s += (3./20.)*temp2;
      temp2*=temp1;
      s += (5./56.)*temp2;
    }
    
    X2(i,0) = s;  
    X2(i,1) = 1.0;
    
    if((using_vertex==true ) && (i == (track.hits.size() - 1)))
    {
      X2(i, 0) /= vertex_sigma_z;
      X2(i, 1) /= vertex_sigma_z;
    }
    else
    {
      X2(i, 0) /= layer_z_resolution[track.hits[i].layer];
      X2(i, 1) /= layer_z_resolution[track.hits[i].layer];
    }
  }
  
  MatrixXf Xt2 = X2.transpose();
  MatrixXf prod2 = Xt2*X2;
  MatrixXf inv2 = prod2.fullPivLu().inverse();
  MatrixXf beta2 = inv2*Xt2*y2;
  
  MatrixXf diff2 = y2 - (X2*beta2);
  MatrixXf chi2_z = (diff2.transpose())*diff2;
  
  float z0 = beta2(1,0);
  float dzdl = beta2(0,0)/sqrt(1. + beta2(0,0)*beta2(0,0));
  
  track.phi = phi;
  track.d = d;
  track.kappa = k;
  track.dzdl = dzdl;
  track.z0 = z0;
  
  if(track.kappa!=0.)
  {
    r=1./track.kappa;
  }
  else
  {
    r=1.0e10;
  }
  
  cx = (track.d+r)*cos(track.phi);
  cy = (track.d+r)*sin(track.phi);
  
  float chi2_tot = 0.;
  for(unsigned int h=0;h<track.hits.size();h++)
  {
    float dx1 = track.hits[h].x - cx;
    float dy1 = track.hits[h].y - cy;
    
    float dx2 = track.hits[h].x + cx;
    float dy2 = track.hits[h].y + cy;
    
    float xydiff1 = sqrt(dx1*dx1 + dy1*dy1) - r;
    float xydiff2 = sqrt(dx2*dx2 + dy2*dy2) - r;
    float xydiff = xydiff2;
    if(fabs(xydiff1) < fabs(xydiff2)){ xydiff = xydiff1; }
    
    float ls_xy = layer_xy_resolution[track.hits[h].layer];
    if((using_vertex == true) && (h == (track.hits.size() - 1)))
    {
      ls_xy = vertex_sigma_xy;
    }
    
    chi2_hit[h] = 0.;
    chi2_hit[h] += xydiff*xydiff/(ls_xy*ls_xy);
    chi2_hit[h] += diff2(h,0)*diff2(h,0);
    
    chi2_tot += chi2_hit[h];
  }
  
  unsigned int deg_of_freedom = 2*track.hits.size() - 5;
  
  if(using_vertex == true)
  {
    track.hits.pop_back();
    chi2_hit.pop_back();
  }
    
  return (chi2_tot)/((double)(deg_of_freedom));
}


