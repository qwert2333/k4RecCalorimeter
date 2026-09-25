#ifndef _FAST_SIPM_DIGI_ALG_C
#define _FAST_SIPM_DIGI_ALG_C

#include "FastSiPMDigiWithEdep.h"

#include <chrono>

// using namespace std;
// using namespace dd4hep;
// using dd4hep::rec::LayeredCalorimeterData;
// using dd4hep::rec::LayeredCalorimeterStruct;

DECLARE_COMPONENT( FastSiPMDigiWithEdep )

FastSiPMDigiWithEdep::FastSiPMDigiWithEdep(const std::string& name, ISvcLocator* svcLoc)
	: Algorithm(name, svcLoc),
    _nEvt(0)
{ 
	// Input collections
    declareProperty("InputSimCaloHitCollection", _inputSimHitCollection, "Handle of the Input SimCaloHit collection");

    // Output collections
    declareProperty("OutputCaloHitCollection", _outputHitCollection, "Handle of Digi CaloHit collection");
    declareProperty("OutputCaloSimLinkCollection", _outputCaloSimLinkCol, "Handle of Sim-CaloHit links");
    declareProperty("OutputCaloMCPLinkCollection", _outputCaloMCPLinkCol, "Handle of Calo-MCParticle links");
}

StatusCode FastSiPMDigiWithEdep::initialize()
{
  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);

  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  // --- Geometry service, cellID decoder
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  m_segmentation = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation();
  if (!m_segmentation) {
    error() << "Readout <<" << m_readoutName << " does not have segmentation" << endmsg;
    return StatusCode::FAILURE;
  }

	// SiPM noise model: Borel distribution
	f_DarkNoise = new TF1("f_DarkNoise", "pow([0]*x, x-1) * exp(-[0]*x) / TMath::Factorial(x)");
	f_DarkNoise->SetParameter(0, _SiPMCT.value());

	// --- Ntuple 
	if(_writeNtuple){
		std::string s_outfile = _filename;
		m_wfile = new TFile(s_outfile.c_str(), "recreate");
		t_DigiHit = new TTree("DigiHit", "DigiHit");
		t_DigiHit->Branch("totE_Truth", &totE_Truth);
		t_DigiHit->Branch("totE_Digi", &totE_Digi);
		t_DigiHit->Branch("digiHit_x", &m_digiHit_x);
		t_DigiHit->Branch("digiHit_y", &m_digiHit_y);
		t_DigiHit->Branch("digiHit_z", &m_digiHit_z);
		t_DigiHit->Branch("digiHit_T", &m_digiHit_T);
		t_DigiHit->Branch("digiHit_E_truth", &m_digiHit_E_truth);
		t_DigiHit->Branch("digiHit_E_digi", &m_digiHit_E_digi);
		t_DigiHit->Branch("digiHit_ADC", &m_digiHit_ADC);
		t_DigiHit->Branch("digiHit_Npe", &m_digiHit_Npe);

		// SiPM
		t_DigiHit->Branch("digiHit_SiPMGain", &m_digiHit_SiPMGain);
		t_DigiHit->Branch("digiHit_SiPMDC", &m_digiHit_SiPMDC);
		t_DigiHit->Branch("digiHit_Pedestal", &m_digiHit_Pedestal);
		t_DigiHit->Branch("digiHit_ADCGain", &m_digiHit_ADCGain);

		// cell ID
		t_DigiHit->Branch("digiHit_stave", &m_digiHit_stave);
		t_DigiHit->Branch("digiHit_sector", &m_digiHit_sector);
		t_DigiHit->Branch("digiHit_rho", &m_digiHit_rho);
		t_DigiHit->Branch("digiHit_theta", &m_digiHit_theta);
		t_DigiHit->Branch("digiHit_phi", &m_digiHit_phi);

	}

	std::cout<<"FastSiPMDigiWithEdep::initialize"<<std::endl;

	return StatusCode::SUCCESS;
}

StatusCode FastSiPMDigiWithEdep::execute(const EventContext&) const
{

	if(_nEvt==0) std::cout<<"FastSiPMDigiWithEdep::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
   	if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

	using Clock = std::chrono::steady_clock;
	using Seconds = std::chrono::duration<double>;
	double loopStepTime = 0.;
	double sortETPairTime = 0.;
	double siPMDigiTime = 0.;

	try{

		const edm4hep::SimCalorimeterHitCollection* simhitCol =  _inputSimHitCollection.get();
		edm4hep::CalorimeterHitCollection* digihitCol = _outputHitCollection.createAndPut();
		auto* caloLinkVec = _outputCaloSimLinkCol.createAndPut();
		auto* caloMCPLinkVec = _outputCaloMCPLinkCol.createAndPut();

    Clear();
		totE_Truth=0;
		totE_Digi=0;
		mean_CT = 0;
		for (int i = 1; i < 10; i++)
		{
			mean_CT += (i-1)*f_DarkNoise->Eval(i);
		}

		// std::cout << "Input sim hit size: " << simhitCol->size() << std::endl;
		
		//Loop in SimHit, digitalize SimHit
		for (const auto& SimHit : *simhitCol) {
			if(!SimHit.isAvailable()) {std::cout<<"Sim hit is not available"<<std::endl; continue;}
			if(SimHit.getEnergy()==0) {std::cout<<"Sim hit energy is 0"<<std::endl; continue;}
      
			double hit_E = SimHit.getEnergy() * dd4hep::GeV;
			totE_Truth += hit_E; 

			// Loop contributions for hit time and Hit-MCParticle connection. 
			std::vector<std::pair<double,double>> step_TimeEnergy;
			MCParticleToEnergyWeightMap MCPEnMap;

			// std::cout << "  Hit energy: " << hit_E << ", step size " << SimHit.contributions_size() << std::endl;

			const auto loopStepStart = Clock::now();
		  for (auto contrib = SimHit.contributions_begin(); contrib != SimHit.contributions_end(); ++contrib) {
      	const double edep = contrib->getEnergy() * dd4hep::GeV;
		  	TVector3 stepPos(contrib->getStepPosition().x * dd4hep::mm, contrib->getStepPosition().y * dd4hep::mm, contrib->getStepPosition().z * dd4hep::mm);


		  	Rndm::Numbers rndm_exp(m_randSvc, Rndm::Exponential(_scintDecaytime.value()));
		  	double init_time = contrib->getTime();
		  	double transport_time = fabs(_outerR.value() - stepPos.Perp()) / (dd4hep::c_light / _refractiveIndex.value());
		  	double scint_time = rndm_exp.shoot();

		  	step_TimeEnergy.emplace_back(init_time + transport_time + scint_time, edep);

		  	auto mcp = contrib->getParticle();
		  	MCPEnMap[mcp] += edep;
		  }
			loopStepTime += Seconds(Clock::now() - loopStepStart).count();

			// std::cout << "  After looping steps: E-T pair size " << step_TimeEnergy.size() << ". Start to sort " << std::endl;

			const auto sortETPairStart = Clock::now();
			std::sort(step_TimeEnergy.begin(), step_TimeEnergy.end(),
          [](auto const &a, auto const &b) {
              return a.first < b.first;
          });
			

			// Time of this hit: first 
			double sumE = 0;
			double hit_T = -1;

			for (auto const &te : step_TimeEnergy) {
				sumE += te.second;
				if (sumE > _Qthfrac*hit_E) {
					hit_T = te.first;
					break;
				}
			}
			sortETPairTime += Seconds(Clock::now() - sortETPairStart).count();

			// std::cout << "  End sort E-T pairs. hit time: " << hit_T << ". Start SiPM digi" << std::endl;

			// #############################################
			// ###########  SiPM Digitization  #############
			// #############################################
			int ScinGen;
			int Npe_SiPM; 
			int integralADC; 
			int darkcounts_CT; 
			double caliE; 
			double calib_const; 
			const auto siPMDigiStart = Clock::now();
			if(_UseDigi){

				Rndm::Numbers rndm_pois(m_randSvc, Rndm::Poisson(hit_E * dd4hep::GeV / dd4hep::MeV  * _CryLY.value()));
				ScinGen = std::round(rndm_pois.shoot());



				// SiPM dark noise and cross talk
				// TODO: the dark count + cross talk and saturation model are not validated yet. 
				//
			  //int darkcounts_mean = Rndm::Numbers(m_randSvc, Rndm::Poisson(_SiPMDCR.value() * _timeWindow.value())).shoot();
			  //darkcounts_CT = 0;
			  //for(int i=0;i<darkcounts_mean;i++)
			  //{
			  //	double darkcounts_rdm = Rndm::Numbers(m_randSvc, Rndm::Flat(0, 1)).shoot();
			  //	int sum_darkcounts = 1;
			  //	if(! (darkcounts_rdm <= f_DarkNoise->Eval(sum_darkcounts)))
			  //	{
			  //		float prob = f_DarkNoise->Eval(sum_darkcounts);
			  //		while(darkcounts_rdm > prob)
			  //		{
			  //			sum_darkcounts++;
			  //			prob += f_DarkNoise->Eval(sum_darkcounts);
			  //		}
			  //	}
			  //	darkcounts_CT += sum_darkcounts;
			  //}

				// std::cout << "Dark count: "<<darkcounts_CT<<std::endl;
				// ScinGen += darkcounts_CT;

				// Npe_SiPM = std::round(_Pixel.value() * (1.0 - TMath::Exp(-ScinGen * 1.0 / _Pixel.value())));
				Npe_SiPM = ScinGen; // Ignore the SiPM saturation effect for now.
				// std::cout<<"After SiPM saturation function: "<<Npe_SiPM<<std::endl;
	

				float ADCMean = Npe_SiPM * _SiPMGainMean.value() * _ADCGainRatio.value() + _Pedestal.value();
				float ADCSigma = std::sqrt(Npe_SiPM * _SiPMGainSigma.value() * _SiPMGainSigma.value() + _PedestalSigma.value() * _PedestalSigma.value());
				integralADC = std::round(Rndm::Numbers(m_randSvc, Rndm::Gauss(ADCMean, ADCSigma)).shoot());

				// std::cout<<"Integrated ADC = "<<integralADC<<" = (Npe "<<Npe_SiPM<<" * Gain "<< _SiPMGainMean.value() * _ADCGainRatio.value() << " + pedestal "<<_Pedestal.value()<<std::endl;

				if(integralADC < 0) integralADC = 0;
				if(integralADC <= _ADCSwitch){  // High gain mode
					calib_const = 1./_SiPMGainMean.value() / _ADCGainRatio.value() / _CryLY.value() * dd4hep::MeV / dd4hep::GeV; 
					caliE = (integralADC-_Pedestal.value()) * calib_const; 
				}
				else if(integralADC > _ADCSwitch && int(integralADC/_ADCGainRatio.value()) <= _ADCSwitch){ // Low gain mode
					ADCMean = Npe_SiPM * _SiPMGainMean.value() + _Pedestal.value();
					ADCSigma = std::sqrt(Npe_SiPM * _SiPMGainSigma.value() * _SiPMGainSigma.value() + _PedestalSigma.value() * _PedestalSigma.value());
					integralADC = std::round(Rndm::Numbers(m_randSvc, Rndm::Gauss(ADCMean, ADCSigma)).shoot());
					if(integralADC < 0) integralADC = 0;

					calib_const = 1./_SiPMGainMean.value() / _CryLY.value()* dd4hep::MeV / dd4hep::GeV; 
					caliE = (integralADC-_Pedestal.value()) * calib_const; 
				}
				else{
					integralADC = _ADC.value()-1;
					calib_const = 1./_SiPMGainMean.value() / _CryLY.value() * dd4hep::MeV / dd4hep::GeV; 
					caliE = (integralADC-_Pedestal.value()) * calib_const; 
				}

				// std::cout << "After ADC gain switch: "<<integralADC<<", calibrated En = " << caliE << " (calib constant " << calib_const << ")" << std::endl;

			}
			else{
				ScinGen = hit_E * dd4hep::GeV / dd4hep::MeV *_CryLY.value();
				Npe_SiPM = ScinGen; 
				integralADC = Npe_SiPM * _SiPMGainMean.value() + _Pedestal.value(); 
				caliE = hit_E; 
			}
			if(Npe_SiPM < _threshold_PE || integralADC<0. || caliE*dd4hep::GeV/dd4hep::MeV < _threshold_PE/_CryLY ) continue; 
			siPMDigiTime += Seconds(Clock::now() - siPMDigiStart).count();

			// std::cout << "  Finish SiPM and ADC digi. Npe = " << Npe_SiPM << ", ADC = " << integralADC << ", energy = " << caliE << std::endl;
			// std::cout << std::endl;
		
		
			// ##################################
			// ####### Some associations  #######
			// ##################################
		
			auto digiHit = digihitCol->create();
			digiHit.setCellID(SimHit.getCellID());
			const double saveVal = (_saveFormat.value() == "Npe") ? Npe_SiPM
			                     : (_saveFormat.value() == "ADC") ? integralADC
			                     : caliE;
			digiHit.setEnergy(saveVal);
			digiHit.setTime(hit_T);
			digiHit.setPosition(SimHit.getPosition());
		
		
			//SimHit - CaloHit association
			auto rel1 = caloLinkVec->create();
			rel1.setFrom(digiHit);
			rel1.setTo(SimHit);
			rel1.setWeight( 1. );		
		
			//MCParticle - CaloHit association
			//float maxMCE = -99.;
			//edm4hep::MCParticle selMCP; 
			for(auto iter : MCPEnMap){
			//if(iter.second>maxMCE){
			//  maxMCE = iter.second;
			//  selMCP = iter.first;
			//}
				auto rel_MCP = caloMCPLinkVec->create();
				rel_MCP.setFrom(digiHit);
				rel_MCP.setTo(iter.first);
				rel_MCP.setWeight(iter.second/SimHit.getEnergy());
			}
		
			totE_Digi += caliE;
		
			if(_writeNtuple){
				m_digiHit_x.push_back(SimHit.getPosition().x);
				m_digiHit_y.push_back(SimHit.getPosition().y);
				m_digiHit_z.push_back(SimHit.getPosition().z);
				m_digiHit_T.push_back(hit_T);
				m_digiHit_E_truth.push_back(hit_E);
				m_digiHit_E_digi.push_back(caliE);
				m_digiHit_ADC.push_back(integralADC);
				m_digiHit_Npe.push_back(Npe_SiPM);
				m_digiHit_SiPMDC.push_back(darkcounts_CT);

				// auto m_decoder = m_segmentation->decoder();
				// auto cellID = SimHit.getCellID();
				// m_digiHit_stave.push_back(static_cast<int>(m_decoder->get(cellID, "stave")));
				// m_digiHit_sector.push_back(static_cast<int>(m_decoder->get(cellID, "sector")));
				// m_digiHit_rho.push_back(static_cast<int>(m_decoder->get(cellID, "rho")));
				// m_digiHit_theta.push_back(static_cast<int>(m_decoder->get(cellID, "theta")));
				// m_digiHit_phi.push_back(static_cast<int>(m_decoder->get(cellID, "phi")));

			}
		}

	}catch(GaudiException &e){
		error()<<"SimCaloHit collection is not available "<<endmsg;
	}
  

	if(_writeNtuple){
		t_DigiHit->Fill();
	}

	debug() << "End Loop: Digitalization!" << endmsg;
	debug() << "Total Truth Energy: " << totE_Truth << endmsg;
	debug() << "Total Digi Energy: " << totE_Digi << endmsg;

	 //const double totalHitCalcTime = loopStepTime + sortETPairTime + siPMDigiTime;
	 //std::cout << "Event " << _nEvt << " timing: "
	 //          << "loop step time = " << loopStepTime * 1000. << " ms, "
	 //          << "sort E-T pair time = " << sortETPairTime * 1000. << " ms, "
	 //          << "SiPM digi time = " << siPMDigiTime * 1000. << " ms, "
	 //          << "total time = " << totalHitCalcTime * 1000. << " ms"
	 //          << std::endl;

	// yyy_enddigi = clock();
	// double duration_digi = double(yyy_enddigi - yyy_start) / CLOCKS_PER_SEC;
	// // 
	// std::ofstream outfile("runtime_ecaldigi.txt", std::ios::app);
	// outfile << _nEvt << "    " << duration_digi << std::endl;
	// outfile.close();

  	_nEvt ++ ;
  	//delete SimHitCol, caloVec, caloLinkVec; 
	return StatusCode::SUCCESS;
}

StatusCode FastSiPMDigiWithEdep::finalize()
{
	if(_writeNtuple && m_wfile){
		m_wfile->cd();
		if(t_DigiHit) t_DigiHit->Write();
		m_wfile->Close();
		delete m_wfile;
		m_wfile = nullptr;
		t_DigiHit = nullptr;
	}
	info() << "Processed " << _nEvt << " events " << endmsg;

	delete f_DarkNoise;
	f_DarkNoise = nullptr;
	return Algorithm::finalize();
}




void FastSiPMDigiWithEdep::Clear() const{
	m_digiHit_x.clear();
	m_digiHit_y.clear();
	m_digiHit_z.clear();
	m_digiHit_T.clear();
	m_digiHit_E_truth.clear();
	m_digiHit_E_digi.clear();
	m_digiHit_ADC.clear();
	m_digiHit_Npe.clear();
	m_digiHit_stave.clear();
	m_digiHit_sector.clear();
	m_digiHit_rho.clear();
	m_digiHit_theta.clear();
	m_digiHit_phi.clear();
	m_digiHit_SiPMGain.clear();
	m_digiHit_SiPMDC.clear();
	m_digiHit_Pedestal.clear();
	m_digiHit_ADCGain.clear();
}

#endif
