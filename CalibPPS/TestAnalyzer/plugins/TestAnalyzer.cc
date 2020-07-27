// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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

class TestAnalyzer : public edm::EDAnalyzer {
public:
  explicit TestAnalyzer(const edm::ParameterSet&);
  ~TestAnalyzer() {};

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
};

TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig) {}


//
// member functions
//

// ------------ method called for each event  ------------
void TestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  ESHandle<PPSAlignmentConfig> config;
  iSetup.get<PPSAlignmentConfigRcd>().get(config);
  std::cout << *config << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestAnalyzer);