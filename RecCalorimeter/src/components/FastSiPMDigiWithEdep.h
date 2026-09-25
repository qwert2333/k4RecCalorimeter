#ifndef _FAST_SIPM_DIGI_ALG_H
#define _FAST_SIPM_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "edm4hep/Vector3f.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableSimCalorimeterHit.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CaloHitMCParticleLinkCollection.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/DetectorData.h>
#include <DD4hep/Segmentations.h> 
#include "k4Interface/IGeoSvc.h"
#include "DD4hep/DD4hepUnits.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"

#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include <cstdlib>

/** @class FastSiPMDigiWithEdep
  *
  *  Simplified algorithm for digitizing the SiPM response using SIM energy deposits and time.
  *  It takes SimCalorimeterHitCollection as input, converting to CalorimeterHitCollection as output.
  *  Options for the output: raw ADC, Number of photon-electrons, or calibrated energy in GeV. 
  *  
  *  Following SiPM effects are included in the digitization:
  *  - Dark count rate (DCR)
  *  - Crosstalk (Xt)
  *  - Saturation (SiPM cell number)

  *  Some ASIC effects are also included, such as:
  *  - Multi-gain mode. 
  *  - Signal-to-noise ratio (SNR)
  *  - Integration time window and threshold

  *  Comparing with SimulateSiPMwithEdep, this simplified version does not simulate the full waveform, 
  *  and ignore all wavelength-dependent effects to fasten the digitization for crystal ECAL with >O(1e4) p.e. per GeV. 

  *  @author Zhiyu Zhao (SJTU) originally developed for CEPC crystal ECAL (now AGORA).
  *  @author Fangyi Guo inherited to FCC software for Grainita ECAL. 
  *  @date   2026-06-22
 */

// #define C 299.79  // unit: mm/ns
// #define PI 3.141592653


class FastSiPMDigiWithEdep : public Gaudi::Algorithm
{
 
public:
    
  FastSiPMDigiWithEdep(const std::string& name, ISvcLocator* svcLoc);
  virtual ~FastSiPMDigiWithEdep() {};

  virtual StatusCode initialize() ;
  virtual StatusCode execute(const EventContext&) const ; 
  virtual StatusCode finalize() ;

private:
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;

  double EnergyDigi(float ScinGen, float sEcalCryMipLY, float sEcalSiPMGainMean, float sEcalSiPMDCR, 
                    TF1* f_SiPMResponse, TF1* f_SiPMSigmaDet, TF1* f_SiPMSigmaRecp, TF1* f_SiPMSigmaRecm, TF1* f_AsymGauss, TF1* f_DarkNoise, TF1* f_ADCNonLin,
                    int& outLO, int& outNDC, int& outNDetPE, float& outPedestal, float& outADC, float& outADCGain);

	void Clear() const;

  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::Segmentation* m_segmentation = nullptr;
  
  // dd4hep::Detector* m_dd4hep;
  // dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	//dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  typedef std::map<const edm4hep::MCParticle, float> MCParticleToEnergyWeightMap;


  TF1* f_DarkNoise = nullptr;
  mutable int mean_CT; 

  // Input names and collections
  Gaudi::Property<std::string> m_readoutName{ this, "ReadOutName", "EcalBarrelCollection" };

  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> _inputSimHitCollection{"GrainitaEcalBarrelRO",
                                                                                 Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> _outputHitCollection{"GrainitaEcalBarrelDigiHit",
                                                                             Gaudi::DataHandle::Writer, this};
  mutable k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> _outputCaloSimLinkCol{
      "GrainitaEcalBarrelDigiHit_SimHit_link", Gaudi::DataHandle::Writer, this};
  mutable k4FWCore::DataHandle<edm4hep::CaloHitMCParticleLinkCollection> _outputCaloMCPLinkCol{
      "GrainitaEcalBarrelDigiHit_MCParticle_link", Gaudi::DataHandle::Writer, this};
  

  mutable Gaudi::Property<int>   _UseDigi{this,   "UseDigi",  1, "If use the digitization model"};

  //Input parameters
  mutable Gaudi::Property<int>   _writeNtuple{this,  "WriteNtuple", 1, "Write ntuple"};
  mutable Gaudi::Property<int>   _Nskip{this,  "SkipEvt", 0, "Skip event"};
  mutable Gaudi::Property<float> _outerR{this,   "OuterRadius", 2131, "ECAL outer radius"};


  // Crystal scintillation parameters
  mutable Gaudi::Property<float> _CryLY{this, 	"LightYield", 10, "Effective light yield (ph/MeV)"};

  // Time parameters
  mutable Gaudi::Property<double> _scintDecaytime{this, "ScintDecaytime", 25000., "scintillation decay time in ns"};
  mutable Gaudi::Property<float> _timeWindow{this, "TimeWindow", 25000., "Time window for digitization in ns"};
  mutable Gaudi::Property<double> _refractiveIndex{this, "RefractiveIndex", 1.6, "Fiber refractive index for transport time estimation"};
  mutable Gaudi::Property<float> _Qthfrac  {this, 	"ChargeThresholdFrac", 0.05, "Charge threshold fraction"};

  // SiPM digitization parameters
  mutable Gaudi::Property<float> _SiPMPDE{this, 	"SiPMPDE", 0.25, "SiPM PDE"};
  mutable Gaudi::Property<float> _SiPMDCR{this, 	"SiPMDCR", 90e-6, "SiPM Dark Count Rate (GHz)"};
  mutable Gaudi::Property<float> _SiPMCT{this, 	"SiPMCT", 0.01, "SiPM crosstalk Probability"};
  mutable Gaudi::Property<int> _Pixel{this, 	"SiPMPixel", 10000, "Number of SiPM pixels"};
  mutable Gaudi::Property<float> _SiPMGainMean{this, 	"SiPMGainMean", 0.0133, "SiPM gain: ADC per p.e. for Gain-1 (ADC)"};
  mutable Gaudi::Property<float> _SiPMGainSigma{this, 	"SiPMGainSigma", 0.001, "Fluctuation of single photoelctron ADC around the mean value for single device (%)"};

  // ADC digitization parameters
  mutable Gaudi::Property<int> _ADC{this, 	"ADC", 4096, "Total ADC conuts for saturation"};
  mutable Gaudi::Property<int> _ADCSwitch{this, 	"ADCSwitch", 4000, "switching point of different gain mode"};
  mutable Gaudi::Property<float> _Pedestal{this, 	"Pedestal", 20, "ADC value of pedestal"};
  mutable Gaudi::Property<float> _ADCGainRatio{this, 	"ADCGainRatio", 30, "Gain ratio"};
  mutable Gaudi::Property<float> _PedestalSigma{this, 	"PedestalResolution", 4, "Pedestal resolution"};


  // Output: threshold and calibration parameters
  mutable Gaudi::Property<float> _threshold_PE{this, 	"ReadoutThreshold", 10., "Threshold for single readout channel. Unit in Npe"};
  mutable Gaudi::Property<std::string> _saveFormat{this, "SaveFormat", "Energy", "Save format in digi hits. Option: Energy, Npe, ADC."};

  // For ntuple

  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};

  typedef std::vector<float> FloatVec;
	mutable int _nEvt ;
	TFile* m_wfile;
	TTree* t_DigiHit;

  mutable double totE_Truth, totE_Digi; 
	mutable FloatVec  m_digiHit_x, m_digiHit_y, m_digiHit_z, 
            m_digiHit_E_truth, m_digiHit_E_digi,
            m_digiHit_T, 
            m_digiHit_ADC, m_digiHit_Npe,
            m_digiHit_stave, m_digiHit_sector, m_digiHit_rho, m_digiHit_theta, m_digiHit_phi;

  mutable FloatVec m_digiHit_SiPMGain, m_digiHit_SiPMDC, m_digiHit_Pedestal, m_digiHit_ADCGain; 

};

#endif
