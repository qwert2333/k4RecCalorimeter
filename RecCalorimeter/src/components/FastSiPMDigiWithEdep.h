#ifndef RECCALORIMETER_FASTSIPMDIGIWITHEDEP_H
#define RECCALORIMETER_FASTSIPMDIGIWITHEDEP_H

#include "GaudiKernel/IRndmGenSvc.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <tuple>

/** Digitize SimCalorimeterHits with a simplified SiPM energy and timing response. */
class FastSiPMDigiWithEdep final
    : public k4FWCore::MultiTransformer<
          std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection,
                     edm4hep::CaloHitMCParticleLinkCollection>(const edm4hep::SimCalorimeterHitCollection&)> {
public:
  FastSiPMDigiWithEdep(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection,
             edm4hep::CaloHitMCParticleLinkCollection>
  operator()(const edm4hep::SimCalorimeterHitCollection& simHits) const override;

private:
  SmartIF<IRndmGenSvc> m_randSvc;
  SmartIF<IGeoSvc> m_geoSvc;

  Gaudi::Property<std::string> m_readoutName{this, "ReadOutName", "EcalBarrelCollection"};
  Gaudi::Property<int> m_useDigi{this, "UseDigi", 1, "Use the digitization model"};
  Gaudi::Property<float> m_outerR{this, "OuterRadius", 2131, "ECAL outer radius in mm"};
  Gaudi::Property<float> m_cryLY{this, "LightYield", 10, "Effective light yield in p.e./MeV"};
  Gaudi::Property<double> m_scintDecaytime{this, "ScintDecaytime", 25000., "Scintillation decay time in ns"};
  Gaudi::Property<float> m_timeWindow{this, "TimeWindow", 25000., "Time window in ns"};
  Gaudi::Property<double> m_refractiveIndex{this, "RefractiveIndex", 1.6};
  Gaudi::Property<float> m_qthfrac{this, "ChargeThresholdFrac", 0.05};

  Gaudi::Property<float> m_sipmPDE{this, "SiPMPDE", 0.25};
  Gaudi::Property<float> m_sipmDCR{this, "SiPMDCR", 90e-6};
  Gaudi::Property<float> m_sipmCT{this, "SiPMCT", 0.01};
  Gaudi::Property<int> m_pixel{this, "SiPMPixel", 10000};
  Gaudi::Property<float> m_sipmGainMean{this, "SiPMGainMean", 0.0133};
  Gaudi::Property<float> m_sipmGainSigma{this, "SiPMGainSigma", 0.001};

  Gaudi::Property<int> m_adc{this, "ADC", 4096};
  Gaudi::Property<int> m_adcSwitch{this, "ADCSwitch", 4000};
  Gaudi::Property<float> m_pedestal{this, "Pedestal", 20};
  Gaudi::Property<float> m_adcGainRatio{this, "ADCGainRatio", 30};
  Gaudi::Property<float> m_pedestalSigma{this, "PedestalResolution", 4};

  Gaudi::Property<float> m_thresholdPE{this, "ReadoutThreshold", 10.};
  Gaudi::Property<std::string> m_saveFormat{this, "SaveFormat", "Energy"};
};

#endif
