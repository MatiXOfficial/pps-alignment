/****************************************************************************
 *
 *  CalibPPS/ESProducers/plugins/PPSAlignmentConfigESSource.cc
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <memory>

//---------------------------------------------------------------------------------------------

class PPSAlignmentConfigESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder
{
public:
    PPSAlignmentConfigESSource(const edm::ParameterSet &);
    ~PPSAlignmentConfigESSource() override = default;

    std::unique_ptr<PPSAlignmentConfig> produce(const PPSAlignmentConfigRcd &);

protected:
	/// sets infinite validity of this data
	void setIntervalFor(const edm::eventsetup::EventSetupRecordKey&,
						const edm::IOVSyncValue&,
						edm::ValidityInterval&) override;

private:
    unsigned int fill;
	unsigned int xangle;
	double beta;
	std::string dataset;

	std::map<unsigned int, std::string> rpTags;

	std::vector<std::string> inputFiles;

	std::map<unsigned int, double> alignmentCorrectionsX, alignmentCorrectionsY;

	bool aligned;

	double n_si;

	SectorConfig sectorConfig45, sectorConfig56;

	std::vector<std::string> matchingReferenceDatasets;
	std::map<unsigned int, SelectionRange> matchingShiftRanges;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges;

	std::map<unsigned int, SelectionRange> alignment_y_ranges;
};

//---------------------------------------------------------------------------------------------

PPSAlignmentConfigESSource::PPSAlignmentConfigESSource(const edm::ParameterSet &iConfig)
{
    rpTags = {
        { 23, "L_2_F" },
		{ 3, "L_1_F" },
		{ 103, "R_1_F" },
		{ 123, "R_2_F" }
    };

    fill = iConfig.getParameter<unsigned int>("fill");
    xangle = iConfig.getParameter<unsigned int>("xangle");
	beta = iConfig.getParameter<double>("beta");
	dataset = iConfig.getParameter<std::string>("dataset");

    const auto &acc = iConfig.getParameter<edm::ParameterSet>("alignment_corrections");
    for (const auto &p : rpTags)
	{
		const auto &ps = acc.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignmentCorrectionsX[p.first] = ps.getParameter<double>("de_x");
		alignmentCorrectionsY[p.first] = ps.getParameter<double>("de_y");
	}

    aligned = iConfig.getParameter<bool>("aligned");

	n_si = iConfig.getParameter<double>("n_si");

	for (std::string sectorName : {"sector_45", "sector_56"})
	{
		const auto &sps = iConfig.getParameter<edm::ParameterSet>(sectorName);
		SectorConfig *sc;
		if (sectorName == "sector_45")
			sc = &sectorConfig45;
		else
			sc = &sectorConfig56;

		sc->cut_h_apply = sps.getParameter<bool>("cut_h_apply");
		sc->cut_h_a = sps.getParameter<double>("cut_h_a");
		sc->cut_h_c = sps.getParameter<double>("cut_h_c");
		sc->cut_h_si = sps.getParameter<double>("cut_h_si");

		sc->cut_v_apply = sps.getParameter<bool>("cut_v_apply");
		sc->cut_v_a = sps.getParameter<double>("cut_v_a");
		sc->cut_v_c = sps.getParameter<double>("cut_v_c");
		sc->cut_v_si = sps.getParameter<double>("cut_v_si");

		sc->nr_x_slice_min = sps.getParameter<double>("nr_x_slice_min");
		sc->nr_x_slice_w = sps.getParameter<double>("nr_x_slice_w");
		sc->nr_x_slice_n = std::ceil((sps.getParameter<double>("nr_x_slice_max") - sc->nr_x_slice_min) / sc->nr_x_slice_w);

		sc->fr_x_slice_min = sps.getParameter<double>("fr_x_slice_min");
		sc->fr_x_slice_w = sps.getParameter<double>("fr_x_slice_w");
		sc->fr_x_slice_n = std::ceil((sps.getParameter<double>("fr_x_slice_max") - sc->fr_x_slice_min) / sc->fr_x_slice_w);
	}

    const auto &c_m = iConfig.getParameter<edm::ParameterSet>("matching");
	matchingReferenceDatasets = c_m.getParameter<std::vector<std::string>>("reference_datasets");

    for (const auto &p : rpTags)
	{
		const auto &ps = c_m.getParameter<edm::ParameterSet>("rp_" + p.second);
		matchingShiftRanges[p.first] = SelectionRange(ps.getParameter<double>("sh_min"), ps.getParameter<double>("sh_max"));
	}

    const auto &c_axo = iConfig.getParameter<edm::ParameterSet>("x_alignment_meth_o");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_axo.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignment_x_meth_o_ranges[p.first] = SelectionRange(ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max"));
	}

    const auto &c_axr = iConfig.getParameter<edm::ParameterSet>("x_alignment_relative");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_axr.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignment_x_relative_ranges[p.first] = SelectionRange(ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max"));
	}

    const auto &c_ay = iConfig.getParameter<edm::ParameterSet>("y_alignment");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_ay.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignment_y_ranges[p.first] = SelectionRange(ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max"));
	}

	setWhatProduced(this);
	findingRecord<PPSAlignmentConfigRcd>();
}

//---------------------------------------------------------------------------------------------

std::unique_ptr<PPSAlignmentConfig> PPSAlignmentConfigESSource::produce(const PPSAlignmentConfigRcd &)
{
    auto p = std::make_unique<PPSAlignmentConfig>();

    p->setFill(fill);
    p->setXangle(xangle);
    p->setBeta(beta);
    p->setDataset(dataset);

    p->setRpTags(rpTags);

    p->setInputFiles(inputFiles);

    p->setAlignmentCorrectionsX(alignmentCorrectionsX);
    p->setAlignmentCorrectionsY(alignmentCorrectionsY);

    p->setAligned(aligned);

    p->setN_si(n_si);

    p->setSectorConfig45(sectorConfig45);
    p->setSectorConfig56(sectorConfig56);

    p->setMatchingReferenceDatasets(matchingReferenceDatasets);
    p->setMatchingShiftRanges(matchingShiftRanges);

    p->setAlignment_x_meth_o_ranges(alignment_x_meth_o_ranges);
    p->setAlignment_x_relative_ranges(alignment_x_relative_ranges);

    p->setAlignment_y_ranges(alignment_y_ranges);

    edm::LogInfo("PPSAlignmentConfigESSource::produce") << "\n" << (*p);

    return p;
}

//---------------------------------------------------------------------------------------------

void PPSAlignmentConfigESSource::setIntervalFor(const edm::eventsetup::EventSetupRecordKey& key,
                                                 const edm::IOVSyncValue& iosv,
                                                 edm::ValidityInterval& oValidity) 
{
	edm::LogInfo("PPSAlignmentConfigESSource")
    	<< ">> PPSAlignmentConfigESSource::setIntervalFor(" << key.name() << ")\n"
    	<< "    run=" << iosv.eventID().run() << ", event=" << iosv.eventID().event();

  	edm::ValidityInterval infinity(iosv.beginOfTime(), iosv.endOfTime());
  	oValidity = infinity;
}

DEFINE_FWK_EVENTSETUP_SOURCE(PPSAlignmentConfigESSource); 