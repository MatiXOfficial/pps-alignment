// -*- C++ -*-
//
// Package:    CalibPPS/TestAnalyzer
// Class:      TestAnalyzer
//
/**\class TestAnalyzer TestAnalyzer.cc CalibPPS/TestAnalyzer/plugins/TestAnalyzer.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mateusz Kocot
//         Created:  Tue, 21 Jul 2020 09:29:28 GMT
//
//

#define THIS_IS_AN_EVENTSETUP_EXAMPLE

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"
//
// class declaration
//

class TestAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TestAnalyzer(const edm::ParameterSet&);
  ~TestAnalyzer() {};

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};


  // ----------member data ---------------------------
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<PPSAlignmentConfig, PPSAlignmentConfigRcd> configToken_;
#endif
};

TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  configToken_ = esConsumes<PPSAlignmentConfig, PPSAlignmentConfigRcd>();
#endif
}

//
// member functions
//

// ------------ method called for each event  ------------
void TestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  auto setup = iSetup.getData(configToken_);
  // std::cout << setup.fill() << std::endl;
#endif
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestAnalyzer);