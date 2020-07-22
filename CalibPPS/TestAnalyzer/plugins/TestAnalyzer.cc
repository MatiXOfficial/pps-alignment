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
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // edm::ESGetToken<PPSAlignmentConfig, PPSAlignmentConfigRcd> configToken_;
};

// TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig) : 
//   configToken_(esConsumes<PPSAlignmentConfig, PPSAlignmentConfigRcd>()) 
//   {}

TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig) {}


//
// member functions
//

// ------------ method called for each event  ------------
void TestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // auto setup = iSetup.getData(configToken_);

  ESHandle<PPSAlignmentConfig> config;
  iSetup.get<PPSAlignmentConfigRcd>().get(config);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestAnalyzer);