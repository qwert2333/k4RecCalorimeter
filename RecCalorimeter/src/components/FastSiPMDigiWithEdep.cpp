#include "FastSiPMDigiWithEdep.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "GaudiKernel/RndmGenerators.h"
#include "k4FWCore/MetadataUtils.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

DECLARE_COMPONENT(FastSiPMDigiWithEdep)

FastSiPMDigiWithEdep::FastSiPMDigiWithEdep(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc, {KeyValue("InputSimCaloHitCollection", "GrainitaEcalBarrelRO")},
                       {KeyValue("OutputCaloHitCollection", "GrainitaEcalBarrelDigiHit"),
                        KeyValue("OutputCaloSimLinkCollection", "GrainitaEcalBarrelDigiHit_SimHit_link"),
                        KeyValue("OutputCaloMCPLinkCollection", "GrainitaEcalBarrelDigiHit_MCParticle_link")}) {}

StatusCode FastSiPMDigiWithEdep::initialize() {
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  const auto& readouts = m_geoSvc->getDetector()->readouts();
  if (readouts.find(m_readoutName.value()) == readouts.end()) {
    error() << "Readout " << m_readoutName << " does not exist" << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_geoSvc->getDetector()->readout(m_readoutName.value()).segmentation().segmentation()) {
    error() << "Readout " << m_readoutName << " does not have segmentation" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_cryLY.value() <= 0 || m_sipmGainMean.value() <= 0 || m_adcGainRatio.value() <= 0 || m_adc.value() <= 0 ||
      m_scintDecaytime.value() <= 0 || m_refractiveIndex.value() <= 0) {
    error() << "LightYield, SiPMGainMean, ADCGainRatio, ADC, ScintDecaytime and RefractiveIndex must be positive"
            << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_saveFormat.value() != "Energy" && m_saveFormat.value() != "Npe" && m_saveFormat.value() != "ADC") {
    error() << "SaveFormat must be Energy, Npe or ADC" << endmsg;
    return StatusCode::FAILURE;
  }

  const auto inputName = inputLocations("InputSimCaloHitCollection")[0];
  const auto outputName = outputLocations("OutputCaloHitCollection")[0];
  if (const auto encoding = k4FWCore::getCellIDEncoding(inputName, this)) {
    k4FWCore::putCellIDEncoding(outputName, *encoding, this);
  }
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection,
           edm4hep::CaloHitMCParticleLinkCollection>
FastSiPMDigiWithEdep::operator()(const edm4hep::SimCalorimeterHitCollection& simHits) const {
  edm4hep::CalorimeterHitCollection digiHits;
  edm4hep::CaloHitSimCaloHitLinkCollection simLinks;
  edm4hep::CaloHitMCParticleLinkCollection mcLinks;

  for (const auto& simHit : simHits) {
    if (!simHit.isAvailable() || simHit.getEnergy() <= 0) {
      continue;
    }

    // EDM4hep energies are in GeV; DD4hep geometry and time calculations use its own units.
    const double hitEnergy = simHit.getEnergy() * dd4hep::GeV;
    std::vector<std::pair<double, double>> timeEnergy;
    std::map<edm4hep::MCParticle, double> mcEnergy;
    for (auto contribution = simHit.contributions_begin(); contribution != simHit.contributions_end(); ++contribution) {
      const double edep = contribution->getEnergy() * dd4hep::GeV;
      const auto& position = contribution->getStepPosition();
      const double radius = std::hypot(position.x, position.y) * dd4hep::mm;
      const double transportTime =
          std::abs(m_outerR.value() * dd4hep::mm - radius) / (dd4hep::c_light / m_refractiveIndex.value());
      Rndm::Numbers scintillation(m_randSvc, Rndm::Exponential(m_scintDecaytime.value()));
      timeEnergy.emplace_back(contribution->getTime() + transportTime + scintillation.shoot(), edep);
      mcEnergy[contribution->getParticle()] += edep;
    }

    std::sort(timeEnergy.begin(), timeEnergy.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    double hitTime = -1;
    double accumulatedEnergy = 0;
    for (const auto& [time, energy] : timeEnergy) {
      accumulatedEnergy += energy;
      if (accumulatedEnergy > m_qthfrac.value() * hitEnergy) {
        hitTime = time;
        break;
      }
    }

    int npe = 0;
    int adc = 0;
    double calibratedEnergy = simHit.getEnergy(); // GeV
    if (m_useDigi.value()) {
      Rndm::Numbers photons(m_randSvc, Rndm::Poisson(hitEnergy / dd4hep::MeV * m_cryLY.value()));
      npe = static_cast<int>(std::round(photons.shoot()));

      const double gainSigma = m_sipmGainSigma.value();
      const double pedestalSigma = m_pedestalSigma.value();
      const double adcSigma = std::sqrt(npe * gainSigma * gainSigma + pedestalSigma * pedestalSigma);
      const double highGainMean = npe * m_sipmGainMean.value() * m_adcGainRatio.value() + m_pedestal.value();
      Rndm::Numbers highGain(m_randSvc, Rndm::Gauss(highGainMean, adcSigma));
      adc = std::max(0, static_cast<int>(std::round(highGain.shoot())));

      double calibration =
          1. / m_sipmGainMean.value() / m_adcGainRatio.value() / m_cryLY.value() * dd4hep::MeV / dd4hep::GeV;
      if (adc > m_adcSwitch.value()) {
        const double lowGainMean = npe * m_sipmGainMean.value() + m_pedestal.value();
        Rndm::Numbers lowGain(m_randSvc, Rndm::Gauss(lowGainMean, adcSigma));
        adc = std::max(0, static_cast<int>(std::round(lowGain.shoot())));
        calibration = 1. / m_sipmGainMean.value() / m_cryLY.value() * dd4hep::MeV / dd4hep::GeV;
      }
      adc = std::min(adc, m_adc.value() - 1);
      calibratedEnergy = (adc - m_pedestal.value()) * calibration;
    } else {
      npe = static_cast<int>(hitEnergy / dd4hep::MeV * m_cryLY.value());
      adc = static_cast<int>(npe * m_sipmGainMean.value() + m_pedestal.value());
    }

    if (npe < m_thresholdPE.value() ||
        calibratedEnergy * dd4hep::GeV / dd4hep::MeV < m_thresholdPE.value() / m_cryLY.value()) {
      continue;
    }

    auto digiHit = digiHits.create();
    digiHit.setCellID(simHit.getCellID());
    digiHit.setPosition(simHit.getPosition());
    digiHit.setTime(hitTime);
    digiHit.setEnergy(m_saveFormat.value() == "Npe" ? npe : m_saveFormat.value() == "ADC" ? adc : calibratedEnergy);

    auto simLink = simLinks.create();
    simLink.setFrom(digiHit);
    simLink.setTo(simHit);
    simLink.setWeight(1.);

    for (const auto& [particle, energy] : mcEnergy) {
      auto mcLink = mcLinks.create();
      mcLink.setFrom(digiHit);
      mcLink.setTo(particle);
      mcLink.setWeight(energy / hitEnergy);
    }
  }

  debug() << "Created " << digiHits.size() << " digitized hits" << endmsg;
  return std::make_tuple(std::move(digiHits), std::move(simLinks), std::move(mcLinks));
}
