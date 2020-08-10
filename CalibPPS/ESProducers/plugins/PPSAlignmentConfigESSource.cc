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
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <memory>

//---------------------------------------------------------------------------------------------

class PPSAlignmentConfigESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder
{
public:
	PPSAlignmentConfigESSource(const edm::ParameterSet &iConfig);
	~PPSAlignmentConfigESSource() override = default;

	std::unique_ptr<PPSAlignmentConfig> produce(const PPSAlignmentConfigRcd &);
	
private:
	void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &key, const edm::IOVSyncValue &iosv, 
	                    edm::ValidityInterval &oValidity) override;

	std::vector<std::string> inputFiles;

	std::map<unsigned int, double> alignmentCorrectionsX, alignmentCorrectionsY;

	bool aligned;

	double n_si;

	SectorConfig sectorConfig45, sectorConfig56;

	std::vector<std::string> matchingReferenceDatasets;
	std::map<unsigned int, SelectionRange> matchingShiftRanges;
	
	std::map<unsigned int, double> yMaxFit;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges;

	std::map<unsigned int, SelectionRange> alignment_y_ranges;

	std::string label;
};

//---------------------------------------------------------------------------------------------

PPSAlignmentConfigESSource::PPSAlignmentConfigESSource(const edm::ParameterSet &iConfig)
{
	for (std::string sectorName : {"sector_45", "sector_56"})
	{
		const auto &sps = iConfig.getParameter<edm::ParameterSet>(sectorName);
		SectorConfig *sc;
		if (sectorName == "sector_45")
			sc = &sectorConfig45;
		else
			sc = &sectorConfig56;

		sc->name = sps.getParameter<std::string>("name");

		for (std::string rpName : {"rp_N", "rp_F"})
		{
			const auto &rps = sps.getParameter<edm::ParameterSet>(rpName);
			RPConfig *rc;
			if (rpName == "rp_N")
				rc = &sc->rp_N;
			else
				rc = &sc->rp_F;

			rc->name = rps.getParameter<std::string>("name");
			rc->id = rps.getParameter<unsigned int>("id");
			rc->position = rps.getParameter<std::string>("position");
			rc->slope = rps.getParameter<double>("slope");
			rc->sh_x = rps.getParameter<double>("sh_x");
		}

		sc->slope = sps.getParameter<double>("slope");
		sc->sh_x_N = sps.getParameter<double>("sh_x_N");

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

	std::map<unsigned int, std::string> rpTags = {
		{ sectorConfig45.rp_F.id, sectorConfig45.rp_F.name },
		{ sectorConfig45.rp_N.id, sectorConfig45.rp_N.name },
		{ sectorConfig56.rp_N.id, sectorConfig56.rp_N.name },
		{ sectorConfig56.rp_F.id, sectorConfig56.rp_F.name }
	};

	const auto &acc = iConfig.getParameter<edm::ParameterSet>("alignment_corrections");
	for (const auto &p : rpTags)
	{
		const auto &ps = acc.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignmentCorrectionsX[p.first] = ps.getParameter<double>("de_x");
		alignmentCorrectionsY[p.first] = ps.getParameter<double>("de_y");
	}

	aligned = iConfig.getParameter<bool>("aligned");

	n_si = iConfig.getParameter<double>("n_si");

	const auto &c_m = iConfig.getParameter<edm::ParameterSet>("matching");
	matchingReferenceDatasets = c_m.getParameter<std::vector<std::string>>("reference_datasets");

	for (const auto &p : rpTags)
	{
		const auto &ps = c_m.getParameter<edm::ParameterSet>("rp_" + p.second);
		matchingShiftRanges[p.first] = {ps.getParameter<double>("sh_min"), ps.getParameter<double>("sh_max")};
	}
	
	const auto &c_mf = iConfig.getParameter<edm::ParameterSet>("y_max_fit");
	for (const auto &p : rpTags)
	{
		yMaxFit[p.first] = c_mf.getParameter<double>("rp_" + p.second);
	}

	const auto &c_axo = iConfig.getParameter<edm::ParameterSet>("x_alignment_meth_o");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_axo.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignment_x_meth_o_ranges[p.first] = {ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max")};
	}

	const auto &c_axr = iConfig.getParameter<edm::ParameterSet>("x_alignment_relative");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_axr.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignment_x_relative_ranges[p.first] = {ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max")};
	}

	const auto &c_ay = iConfig.getParameter<edm::ParameterSet>("y_alignment");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_ay.getParameter<edm::ParameterSet>("rp_" + p.second);
		alignment_y_ranges[p.first] = {ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max")};
	}

	label = iConfig.getParameter<std::string>("label");
	setWhatProduced(this, label);
	findingRecord<PPSAlignmentConfigRcd>();
}

//---------------------------------------------------------------------------------------------

std::unique_ptr<PPSAlignmentConfig> PPSAlignmentConfigESSource::produce(const PPSAlignmentConfigRcd &)
{
	auto p = std::make_unique<PPSAlignmentConfig>();

	p->setSectorConfig45(sectorConfig45);
	p->setSectorConfig56(sectorConfig56);

	p->setInputFiles(inputFiles);

	p->setAlignmentCorrectionsX(alignmentCorrectionsX);
	p->setAlignmentCorrectionsY(alignmentCorrectionsY);

	p->setAligned(aligned);

	p->setN_si(n_si);

	p->setMatchingReferenceDatasets(matchingReferenceDatasets);
	p->setMatchingShiftRanges(matchingShiftRanges);
	
	p->setYMaxFit(yMaxFit);

	p->setAlignment_x_meth_o_ranges(alignment_x_meth_o_ranges);
	p->setAlignment_x_relative_ranges(alignment_x_relative_ranges);

	p->setAlignment_y_ranges(alignment_y_ranges);

	edm::LogInfo("produce_" + label) << "\n" << (*p);

	return p;
}

//---------------------------------------------------------------------------------------------

void PPSAlignmentConfigESSource::setIntervalFor(const edm::eventsetup::EventSetupRecordKey& key,
												const edm::IOVSyncValue& iosv,
												edm::ValidityInterval& oValidity) 
{
	edm::LogInfo("PPSAlignmentConfigESSource")
	<< ">> PPSAlignmentConfigESSource_setIntervalFor(" << key.name() << ")\n"
	<< "    run=" << iosv.eventID().run() << ", event=" << iosv.eventID().event();

	edm::ValidityInterval infinity(iosv.beginOfTime(), iosv.endOfTime());
	oValidity = infinity;
}

DEFINE_FWK_EVENTSETUP_SOURCE(PPSAlignmentConfigESSource); 