//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: Henri Payno and Yann Perrot
//
//
/// \file ClusteringAlgorithm.hh
/// \brief Definition of the ClusteringAlgorithm class

#ifndef ClusteringAlgorithm_H
#define ClusteringAlgorithm_H 1

#include "ClusterSBPoints_.hh"
#include "SBPoint_.hh"

#include <map>

class ClusteringAlgorithm
{
public:

  ClusteringAlgorithm(int pEps, int pMinPts,
      double pEMinDamage, double pEMaxDamage);
  ~ClusteringAlgorithm();

  // Get Set methods
  int GetEps()
  {
    return fEps;
  };
  void SetEps(int val)
  {
    fEps=val;
  };
  int GetMinPts()
  {
    return fMinPts;
  };
  void SetMinPts(int val)
  {
    fMinPts=val;
  };
  double GetSPointsProb()
  {
    return fSPointsProb;
  };
  void SetSPointsProb(double val)
  {
    fSPointsProb=val;
  };
  double GetEMinDamage()
  {
    return fEMinDamage;
  };
  void SetEMinDamage(double val)
  {
    fEMinDamage=val;
  };
  double GetEMaxDamage()
  {
    return fEMaxDamage;
  };
  void SetEMaxDamage(double val)
  {
    fEMaxDamage=val;
  };

  // Register a damage (position, edep)
  void RegisterDamage(int, double, int);

  // Clustering Algorithm
  std::map<int,int> RunClustering();

  // Clean all data structures
  void  Purge();

  // Return the number of simple break
  int GetSSB() const;
  // Return the number of complex simple break
  int GetComplexSSB() const;
  // Return the number of double strand break
  int GetDSB() const;
  // Return a map representing cluster size distribution
  // first int : cluster size (1 = SSB)
  // second int : counts

  int GetTotalSB()
  {
    return fpSetOfPoints.size();
  }
  std::map<int,int> GetClusterSizeDistribution();

  bool hasRegistered {false};
  // Functions to check if SB candidate
  bool IsEdepSufficient(double);


private:

  // Check if a SB point can be merged to a cluster, and do it
  bool FindCluster(SBPoint* pPt);
  // Check if two points can be merged
  bool AreOnTheSameCluster(int, int,int);
  // Merge clusters
  void MergeClusters();
  // Add SSB to clusters
  void IncludeUnassociatedPoints();

  // Parameters to run clustering algorithm
  int fEps;         // distance to merge SBPoints
  int fMinPts;         // number of SBPoints to create a cluster
  double fSPointsProb; // probability for a point to be in the sensitive area
  double fEMinDamage;  // min energy to create a damage
  double fEMaxDamage;  // energy to have a probability to create a damage = 1

  // Data structure containing all SB points
  std::vector<SBPoint*> fpSetOfPoints;
  // Data structure containing all clusters
  std::vector<ClusterSBPoints*> fpClusters;
  // ID of the next SB point
  unsigned int fNextSBPointID;


};

#endif

