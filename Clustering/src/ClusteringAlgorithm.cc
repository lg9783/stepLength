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
/// \file ClusteringAlgo.cc
/// \brief Implementation of the ClustreringAlgo class

#include "ClusteringAlgorithm.hh"

// #include "TRandom.h"

#include <map>
#include <stdlib.h> 

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgorithm::ClusteringAlgorithm(int pEps,int pMinPts,double pEMinDamage, double pEMaxDamage)
:fEps(pEps),fMinPts(pMinPts),fEMinDamage(pEMinDamage),fEMaxDamage(pEMaxDamage)
{
  fNextSBPointID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgorithm::~ClusteringAlgorithm()
{
  Purge();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Random sampling in energy
bool ClusteringAlgorithm::IsEdepSufficient(double pEdep)
{
  // if(pEdep<fEMinDamage)
  // {
  //   return false;
  // }

  // else if(pEdep>fEMaxDamage)
  // {
  //   return true;
  // }
  // else
  // {
  //   double proba = (pEdep - fEMinDamage)/
  //       (fEMaxDamage-fEMinDamage);
  //   return (proba>R.Uniform());
    
  // }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add an event interaction to the unregistered damage if
// good conditions (pos and energy) are met
//

void ClusteringAlgorithm::RegisterDamage(int pPos,double pEdep, int strand)
{
      fpSetOfPoints.push_back( new SBPoint(fNextSBPointID++, pPos, pEdep, strand));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

map<int,int> ClusteringAlgorithm::RunClustering()
{
  // quick sort style
  // create cluster
  std::vector<SBPoint*>::iterator itVisitorPt, itObservedPt;
  for(itVisitorPt = fpSetOfPoints.begin();
      itVisitorPt != fpSetOfPoints.end();
      ++itVisitorPt  )
  {
    itObservedPt = itVisitorPt;
    itObservedPt ++;

    while(itObservedPt != fpSetOfPoints.end() )
    {
      // if at least one of the two points has not a cluster
      if(!((*itObservedPt)->HasCluster() && (*itVisitorPt)->HasCluster()))
      {
        if(AreOnTheSameCluster( (*itObservedPt)->GetPosition(),
                                (*itVisitorPt)->GetPosition(),fEps))
        {
          // if none has a cluster. Create a new one
          if(!(*itObservedPt)->HasCluster() && !(*itVisitorPt)->HasCluster())
          {
            // create the new cluster
            set<SBPoint*> clusterPoints;
            clusterPoints.insert((*itObservedPt));
            clusterPoints.insert((*itVisitorPt));
            ClusterSBPoints* lCluster = new ClusterSBPoints(clusterPoints);
            assert(lCluster);
            fpClusters.push_back(lCluster);
            assert(lCluster);
            // inform SB point that they are part of a cluster now
            assert(lCluster);
            (*itObservedPt)->SetCluster(lCluster);
            assert(lCluster);
            (*itVisitorPt)->SetCluster(lCluster);
          }else
          {
            // add the point to the existing cluster
            if((*itObservedPt)->HasCluster())
            {
              (*itObservedPt)->GetCluster()->AddSBPoint((*itVisitorPt));
              (*itVisitorPt)->SetCluster((*itObservedPt)->GetCluster());
            }

            if((*itVisitorPt)->HasCluster())
            {
              (*itVisitorPt)->GetCluster()->AddSBPoint((*itObservedPt));
              (*itObservedPt)->SetCluster((*itVisitorPt)->GetCluster());
            }
          }
        }
      }
      ++itObservedPt;
    }
  }

  // associate isolated points and merge clusters
  IncludeUnassociatedPoints();
  MergeClusters();

  // return cluster size distribution
  return GetClusterSizeDistribution();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// try to merge cluster between them, based on the distance between barycenters
void ClusteringAlgorithm::MergeClusters()
{
  std::vector<ClusterSBPoints*>::iterator itCluster1, itCluster2;
  for(itCluster1 = fpClusters.begin();
      itCluster1 != fpClusters.end();
      ++itCluster1)
  {
    int baryCenterClust1 = (*itCluster1)->GetBarycenter();
    itCluster2 = itCluster1;
    itCluster2++;
    while(itCluster2 != fpClusters.end())
    {
      int baryCenterClust2 = (*itCluster2)->GetBarycenter();
      // if we can merge both cluster
      if(AreOnTheSameCluster(baryCenterClust1, baryCenterClust2,fEps))
      {
        (*itCluster1)->MergeWith(*itCluster2);
        delete *itCluster2;
        fpClusters.erase(itCluster2);
        return MergeClusters();
      }else
      {
        itCluster2++;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgorithm::IncludeUnassociatedPoints()
{
  std::vector<SBPoint*>::iterator itVisitorPt;
  int nbPtSansCluster = 0;
  // Associate all point not in a cluster if possible ( to the first found cluster)
  for(itVisitorPt = fpSetOfPoints.begin();
      itVisitorPt != fpSetOfPoints.end();
      ++itVisitorPt)
  {
    if(!(*itVisitorPt)->HasCluster())
    {
      nbPtSansCluster ++;
      FindCluster(*itVisitorPt);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusteringAlgorithm::FindCluster(SBPoint* pPt)
{
  assert(!pPt->HasCluster());
  std::vector<ClusterSBPoints*>::iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    //if((*itCluster)->hasIn(pPt, fEps))
    if((*itCluster)->HasInBarycenter(pPt, fEps))
    {
      (*itCluster)->AddSBPoint(pPt);
      return true;
    }
  }  
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusteringAlgorithm::AreOnTheSameCluster(int copy1,
    int copy2, int pMinBP)
{

  if(abs(copy1-copy2)<= pMinBP)
  {
    return true;
  }else
  {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ClusteringAlgorithm::GetSSB() const
{
  int nbSSB = 0;
  std::vector<SBPoint*>::const_iterator itSDSPt;
  for(itSDSPt = fpSetOfPoints.begin();
      itSDSPt != fpSetOfPoints.end();
      ++itSDSPt)
  {
    if(!(*itSDSPt)->HasCluster())
    {
      nbSSB++;
    }
  }
  return nbSSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ClusteringAlgorithm::GetComplexSSB() const
{
  int nbSSB = 0;
  std::vector<ClusterSBPoints*>::const_iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    if((*itCluster)->IsSSB())
    {
      nbSSB ++;
    }
  }
  return nbSSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ClusteringAlgorithm::GetDSB() const
{
  int nbDSB = 0;
  std::vector<ClusterSBPoints*>::const_iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    if((*itCluster)->IsDSB())
    {
      nbDSB ++;
    }
  }
  return nbDSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

map<int,int> ClusteringAlgorithm::GetClusterSizeDistribution()
{
  std::map<int,int>  sizeDistribution;
  sizeDistribution[1]=GetSSB();
  std::vector<ClusterSBPoints*>::const_iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      itCluster++)
  {
    sizeDistribution[(*itCluster)->GetSize()]++;
  }
  return sizeDistribution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgorithm::Purge()
{
  fNextSBPointID = 0;
  std::vector<ClusterSBPoints*>::iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    delete *itCluster;
    *itCluster = NULL;
  }  
  fpClusters.clear();
  std::vector<SBPoint*>::iterator itPt;
  for(itPt = fpSetOfPoints.begin();
      itPt != fpSetOfPoints.end();
      ++itPt)
  {
    delete *itPt;
    *itPt = NULL;
  }
  fpSetOfPoints.clear();
}

