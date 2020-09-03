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
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <string>
#include <vector>
#include <map>
#include <memory>

//---------------------------------------------------------------------------------------------

class PPSAlignmentConfigESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder
{
public:
	PPSAlignmentConfigESSource(const edm::ParameterSet &iConfig);
	~PPSAlignmentConfigESSource() override = default;

	std::unique_ptr<PPSAlignmentConfig> produce(const PPSAlignmentConfigRcd &);
	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
	
private:
	void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &key, const edm::IOVSyncValue &iosv, 
	                    edm::ValidityInterval &oValidity) override;

	std::vector<std::string> sequence;

	SectorConfig sectorConfig45, sectorConfig56;

	std::map<unsigned int, double> alignmentCorrectionsX, alignmentCorrectionsY;

	double x_ali_sh_step;

	double y_mode_sys_unc;
	double chiSqThreshold;
	double y_mode_unc_max_valid;
	double y_mode_max_valid;

	double n_si;

	std::vector<std::string> matchingReferenceDatasets;
	std::map<unsigned int, SelectionRange> matchingShiftRanges;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges;

	std::map<unsigned int, SelectionRange> alignment_y_ranges;

	Binning binning;

	std::string label;
};

//---------------------------------------------------------------------------------------------

PPSAlignmentConfigESSource::PPSAlignmentConfigESSource(const edm::ParameterSet &iConfig)
{
	sequence = iConfig.getParameter<std::vector<std::string>>("sequence");

	sectorConfig45.name = "sector 45";

	sectorConfig45.rp_N.position = "N";
	sectorConfig45.rp_F.position = "F";

	sectorConfig56.name = "sector 56";

	sectorConfig56.rp_N.position = "N";
	sectorConfig56.rp_F.position = "F";

	for (std::string sectorName : {"sector_45", "sector_56"})
	{
		const auto &sps = iConfig.getParameter<edm::ParameterSet>(sectorName);
		SectorConfig *sc;
		if (sectorName == "sector_45")
			sc = &sectorConfig45;
		else
			sc = &sectorConfig56;

		const auto &rnps = sps.getParameter<edm::ParameterSet>("rp_N");
		sc->rp_N.name = rnps.getParameter<std::string>("name");
		sc->rp_N.id = rnps.getParameter<int>("id");

		sc->rp_N.slope = rnps.getParameter<double>("slope");
		sc->rp_N.sh_x = rnps.getParameter<double>("sh_x");

		sc->rp_N.x_min_mode = rnps.getParameter<double>("x_min_mode");
		sc->rp_N.x_max_mode = rnps.getParameter<double>("x_max_mode");
		sc->rp_N.y_cen_add = rnps.getParameter<double>("y_cen_add");
		sc->rp_N.y_width_mult = rnps.getParameter<double>("y_width_mult");

		const auto &rfps = sps.getParameter<edm::ParameterSet>("rp_F");
		sc->rp_F.name = rfps.getParameter<std::string>("name");
		sc->rp_F.id = rfps.getParameter<int>("id");

		sc->rp_F.slope = rfps.getParameter<double>("slope");
		sc->rp_F.sh_x = rfps.getParameter<double>("sh_x");

		sc->rp_F.x_min_mode = rfps.getParameter<double>("x_min_mode");
		sc->rp_F.x_max_mode = rfps.getParameter<double>("x_max_mode");
		sc->rp_F.y_cen_add = rfps.getParameter<double>("y_cen_add");
		sc->rp_F.y_width_mult = rfps.getParameter<double>("y_width_mult");

		sc->slope = sps.getParameter<double>("slope");

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
		{ sectorConfig45.rp_F.id, "rp_L_F" },
		{ sectorConfig45.rp_N.id, "rp_L_N" },
		{ sectorConfig56.rp_N.id, "rp_R_N" },
		{ sectorConfig56.rp_F.id, "rp_R_F" }
	};

	const auto &acc = iConfig.getParameter<edm::ParameterSet>("alignment_corrections");
	for (const auto &p : rpTags)
	{
		const auto &ps = acc.getParameter<edm::ParameterSet>(p.second);
		alignmentCorrectionsX[p.first] = ps.getParameter<double>("de_x");
		alignmentCorrectionsY[p.first] = ps.getParameter<double>("de_y");
	}

	x_ali_sh_step = iConfig.getParameter<double>("x_ali_sh_step");

	y_mode_sys_unc = iConfig.getParameter<double>("y_mode_sys_unc");
	chiSqThreshold = iConfig.getParameter<double>("chiSqThreshold");
	y_mode_unc_max_valid = iConfig.getParameter<double>("y_mode_unc_max_valid");
	y_mode_max_valid = iConfig.getParameter<double>("y_mode_max_valid");

	n_si = iConfig.getParameter<double>("n_si");

	const auto &c_m = iConfig.getParameter<edm::ParameterSet>("matching");
	matchingReferenceDatasets = c_m.getParameter<std::vector<std::string>>("reference_datasets");

	for (const auto &p : rpTags)
	{
		const auto &ps = c_m.getParameter<edm::ParameterSet>(p.second);
		matchingShiftRanges[p.first] = {ps.getParameter<double>("sh_min"), ps.getParameter<double>("sh_max")};
	}

	const auto &c_axo = iConfig.getParameter<edm::ParameterSet>("x_alignment_meth_o");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_axo.getParameter<edm::ParameterSet>(p.second);
		alignment_x_meth_o_ranges[p.first] = {ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max")};
	}

	const auto &c_axr = iConfig.getParameter<edm::ParameterSet>("x_alignment_relative");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_axr.getParameter<edm::ParameterSet>(p.second);
		alignment_x_relative_ranges[p.first] = {ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max")};
	}

	const auto &c_ay = iConfig.getParameter<edm::ParameterSet>("y_alignment");
	for (const auto &p : rpTags)
	{
		const auto &ps = c_ay.getParameter<edm::ParameterSet>(p.second);
		alignment_y_ranges[p.first] = {ps.getParameter<double>("x_min"), ps.getParameter<double>("x_max")};
	}

	const auto &bps = iConfig.getParameter<edm::ParameterSet>("binning");
	binning.bin_size_x = bps.getParameter<double>("bin_size_x");
	binning.n_bins_x = bps.getParameter<unsigned int>("n_bins_x");
	binning.pixel_x_offset = bps.getParameter<double>("pixel_x_offset");
	binning.n_bins_y = bps.getParameter<unsigned int>("n_bins_y");
	binning.y_min = bps.getParameter<double>("y_min");
	binning.y_max = bps.getParameter<double>("y_max");

	label = iConfig.getParameter<std::string>("label");
	setWhatProduced(this, label);
	findingRecord<PPSAlignmentConfigRcd>();
}

//---------------------------------------------------------------------------------------------

std::unique_ptr<PPSAlignmentConfig> PPSAlignmentConfigESSource::produce(const PPSAlignmentConfigRcd &)
{
	auto p = std::make_unique<PPSAlignmentConfig>();

	p->setSequence(sequence);

	p->setSectorConfig45(sectorConfig45);
	p->setSectorConfig56(sectorConfig56);

	p->setAlignmentCorrectionsX(alignmentCorrectionsX);
	p->setAlignmentCorrectionsY(alignmentCorrectionsY);

	p->setX_ali_sh_step(x_ali_sh_step);

	p->setY_mode_sys_unc(y_mode_sys_unc);
	p->setChiSqThreshold(chiSqThreshold);
	p->setY_mode_unc_max_valid(y_mode_unc_max_valid);
	p->setY_mode_max_valid(y_mode_max_valid);

	p->setN_si(n_si);

	p->setMatchingReferenceDatasets(matchingReferenceDatasets);
	p->setMatchingShiftRanges(matchingShiftRanges);

	p->setAlignment_x_meth_o_ranges(alignment_x_meth_o_ranges);
	p->setAlignment_x_relative_ranges(alignment_x_relative_ranges);

	p->setAlignment_y_ranges(alignment_y_ranges);

	p->setBinning(binning);

	edm::LogInfo("produce") << "\n" << (label.empty() ? "empty label" : "label = " + label) << ":\n\n" << (*p);

	return p;
}

//---------------------------------------------------------------------------------------------

/**********************
defaults for:
	* fill = 7334
	* xangle = 160
	* beta = 0.30
	* data not already aligned
	* skeleton from 2018 config_base.py
**********************/
void PPSAlignmentConfigESSource::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;

	desc.add<std::string>("label", "");

	desc.add<std::vector<std::string>>("sequence", {});

	// sector_45
	{
		edm::ParameterSetDescription sector45;

		edm::ParameterSetDescription rp_N;
		rp_N.add<std::string>("name", "L_1_F");
		rp_N.add<int>("id", 3);

		rp_N.add<double>("slope", 0.19);
		rp_N.add<double>("sh_x", -3.6);

		rp_N.add<double>("x_min_mode", 2.);
		rp_N.add<double>("x_max_mode", 7.0);
		rp_N.add<double>("y_cen_add", -0.3);
		rp_N.add<double>("y_width_mult", 1.1);
		sector45.add<edm::ParameterSetDescription>("rp_N", rp_N);

		edm::ParameterSetDescription rp_F;
		rp_F.add<std::string>("name", "L_2_F");
		rp_F.add<int>("id", 23);
		
		rp_F.add<double>("slope", 0.19);
		rp_F.add<double>("sh_x", -42.);

		rp_F.add<double>("x_min_mode", 2.);
		rp_F.add<double>("x_max_mode", 7.5);
		rp_F.add<double>("y_cen_add", -0.3);
		rp_F.add<double>("y_width_mult", 1.1);
		sector45.add<edm::ParameterSetDescription>("rp_F", rp_F);

		sector45.add<double>("slope", 0.006);
		sector45.add<bool>("cut_h_apply", true);
		sector45.add<double>("cut_h_a", -1.);
		sector45.add<double>("cut_h_c", -38.55);
		sector45.add<double>("cut_h_si", 0.2);
		sector45.add<bool>("cut_v_apply", true);
		sector45.add<double>("cut_v_a", -1.07);
		sector45.add<double>("cut_v_c", 1.63);
		sector45.add<double>("cut_v_si", 0.15);
		sector45.add<double>("nr_x_slice_min", 7.);
		sector45.add<double>("nr_x_slice_max", 19.);
		sector45.add<double>("nr_x_slice_w", 0.2);
		sector45.add<double>("fr_x_slice_min", 46.);
		sector45.add<double>("fr_x_slice_max", 58.);
		sector45.add<double>("fr_x_slice_w", 0.2);

		desc.add<edm::ParameterSetDescription>("sector_45", sector45);
	}

	// sector_56
	{
		edm::ParameterSetDescription sector56;

		edm::ParameterSetDescription rp_N;
		rp_N.add<std::string>("name", "R_1_F");
		rp_N.add<int>("id", 103);
		
		rp_N.add<double>("slope", 0.40);
		rp_N.add<double>("sh_x", -2.8);

		rp_N.add<double>("x_min_mode", 2.);
		rp_N.add<double>("x_max_mode", 7.4);
		rp_N.add<double>("y_cen_add", -0.8);
		rp_N.add<double>("y_width_mult", 1.0);
		sector56.add<edm::ParameterSetDescription>("rp_N", rp_N);

		edm::ParameterSetDescription rp_F;
		rp_F.add<std::string>("name", "R_2_F");
		rp_F.add<int>("id", 123);
		
		rp_F.add<double>("slope", 0.39);
		rp_F.add<double>("sh_x", -41.9);

		rp_F.add<double>("x_min_mode", 2.);
		rp_F.add<double>("x_max_mode", 8.0);
		rp_F.add<double>("y_cen_add", -0.8);
		rp_F.add<double>("y_width_mult", 1.0);
		sector56.add<edm::ParameterSetDescription>("rp_F", rp_F);

		sector56.add<double>("slope", -0.015);
		sector56.add<bool>("cut_h_apply", true);
		sector56.add<double>("cut_h_a", -1.);
		sector56.add<double>("cut_h_c", -39.26);
		sector56.add<double>("cut_h_si", 0.2);
		sector56.add<bool>("cut_v_apply", true);
		sector56.add<double>("cut_v_a", -1.07);
		sector56.add<double>("cut_v_c", 1.49);
		sector56.add<double>("cut_v_si", 0.15);
		sector56.add<double>("nr_x_slice_min", 6.);
		sector56.add<double>("nr_x_slice_max", 17.);
		sector56.add<double>("nr_x_slice_w", 0.2);
		sector56.add<double>("fr_x_slice_min", 45.);
		sector56.add<double>("fr_x_slice_max", 57.);
		sector56.add<double>("fr_x_slice_w", 0.2);

		desc.add<edm::ParameterSetDescription>("sector_56", sector56);
	}

	// alignment_corrections
	{
		edm::ParameterSetDescription alignmentCorrections;

		edm::ParameterSetDescription rpLF;
		rpLF.add<double>("de_x", 0.);
		rpLF.add<double>("de_y", 0.);
		alignmentCorrections.add<edm::ParameterSetDescription>("rp_L_F", rpLF);

		edm::ParameterSetDescription rpLN;
		rpLN.add<double>("de_x", 0.);
		rpLN.add<double>("de_y", 0.);
		alignmentCorrections.add<edm::ParameterSetDescription>("rp_L_N", rpLN);

		edm::ParameterSetDescription rpRN;
		rpRN.add<double>("de_x", 0.);
		rpRN.add<double>("de_y", 0.);
		alignmentCorrections.add<edm::ParameterSetDescription>("rp_R_N", rpRN);

		edm::ParameterSetDescription rpRF;
		rpRF.add<double>("de_x", 0.);
		rpRF.add<double>("de_y", 0.);
		alignmentCorrections.add<edm::ParameterSetDescription>("rp_R_F", rpRF);

		desc.add<edm::ParameterSetDescription>("alignment_corrections", alignmentCorrections);
	}

	desc.add<double>("x_ali_sh_step", 0.01);

	desc.add<double>("y_mode_sys_unc", 0.03);
	desc.add<double>("chiSqThreshold", 50.);
	desc.add<double>("y_mode_unc_max_valid", 5.);
	desc.add<double>("y_mode_max_valid", 20.);

	desc.add<double>("n_si", 4.);

	// matching
	{
		edm::ParameterSetDescription matching;
		std::vector<std::string> referenceDatasets;
		matching.add<std::vector<std::string>>("reference_datasets", referenceDatasets);

		edm::ParameterSetDescription rpLF;
		rpLF.add<double>("sh_min", -43.);
		rpLF.add<double>("sh_max", -41.);
		matching.add<edm::ParameterSetDescription>("rp_L_F", rpLF);

		edm::ParameterSetDescription rpLN;
		rpLN.add<double>("sh_min", -4.2);
		rpLN.add<double>("sh_max", -2.4);
		matching.add<edm::ParameterSetDescription>("rp_L_N", rpLN);

		edm::ParameterSetDescription rpRN;
		rpRN.add<double>("sh_min", -3.6);
		rpRN.add<double>("sh_max", -1.8);
		matching.add<edm::ParameterSetDescription>("rp_R_N", rpRN);

		edm::ParameterSetDescription rpRF;
		rpRF.add<double>("sh_min", -43.2);
		rpRF.add<double>("sh_max", -41.2);
		matching.add<edm::ParameterSetDescription>("rp_R_F", rpRF);

		desc.add<edm::ParameterSetDescription>("matching", matching);
	}

	// x alignment meth o
	{
		edm::ParameterSetDescription x_alignment_meth_o;

		edm::ParameterSetDescription rpLF;
		rpLF.add<double>("x_min", 47.);
		rpLF.add<double>("x_max", 56.5);
		x_alignment_meth_o.add<edm::ParameterSetDescription>("rp_L_F", rpLF);

		edm::ParameterSetDescription rpLN;
		rpLN.add<double>("x_min", 9.);
		rpLN.add<double>("x_max", 18.5);
		x_alignment_meth_o.add<edm::ParameterSetDescription>("rp_L_N", rpLN);

		edm::ParameterSetDescription rpRN;
		rpRN.add<double>("x_min", 7.);
		rpRN.add<double>("x_max", 15.);
		x_alignment_meth_o.add<edm::ParameterSetDescription>("rp_R_N", rpRN);

		edm::ParameterSetDescription rpRF;
		rpRF.add<double>("x_min", 46.);
		rpRF.add<double>("x_max", 54.);
		x_alignment_meth_o.add<edm::ParameterSetDescription>("rp_R_F", rpRF);

		desc.add<edm::ParameterSetDescription>("x_alignment_meth_o", x_alignment_meth_o);
	}

	// x alignment relative
	{
		edm::ParameterSetDescription x_alignment_relative;

		edm::ParameterSetDescription rpLF;
		rpLF.add<double>("x_min", 0.);
		rpLF.add<double>("x_max", 0.);
		x_alignment_relative.add<edm::ParameterSetDescription>("rp_L_F", rpLF);

		edm::ParameterSetDescription rpLN;
		rpLN.add<double>("x_min", 7.5);
		rpLN.add<double>("x_max", 12.);
		x_alignment_relative.add<edm::ParameterSetDescription>("rp_L_N", rpLN);

		edm::ParameterSetDescription rpRN;
		rpRN.add<double>("x_min", 6.);
		rpRN.add<double>("x_max", 10.);
		x_alignment_relative.add<edm::ParameterSetDescription>("rp_R_N", rpRN);

		edm::ParameterSetDescription rpRF;
		rpRF.add<double>("x_min", 0.);
		rpRF.add<double>("x_max", 0.);
		x_alignment_relative.add<edm::ParameterSetDescription>("rp_R_F", rpRF);

		desc.add<edm::ParameterSetDescription>("x_alignment_relative", x_alignment_relative);
	}

	// y alignment
	{
		edm::ParameterSetDescription y_alignment;

		edm::ParameterSetDescription rpLF;
		rpLF.add<double>("x_min", 44.5);
		rpLF.add<double>("x_max", 49.);
		y_alignment.add<edm::ParameterSetDescription>("rp_L_F", rpLF);

		edm::ParameterSetDescription rpLN;
		rpLN.add<double>("x_min", 6.7);
		rpLN.add<double>("x_max", 11.);
		y_alignment.add<edm::ParameterSetDescription>("rp_L_N", rpLN);

		edm::ParameterSetDescription rpRN;
		rpRN.add<double>("x_min", 5.9);
		rpRN.add<double>("x_max", 10.);
		y_alignment.add<edm::ParameterSetDescription>("rp_R_N", rpRN);

		edm::ParameterSetDescription rpRF;
		rpRF.add<double>("x_min", 44.5);
		rpRF.add<double>("x_max", 49.);
		y_alignment.add<edm::ParameterSetDescription>("rp_R_F", rpRF);

		desc.add<edm::ParameterSetDescription>("y_alignment", y_alignment);
	}

	// binning
	{
		edm::ParameterSetDescription binning;

		binning.add<double>("bin_size_x", 142.3314E-3);
		binning.add<unsigned int>("n_bins_x", 210);
		binning.add<double>("pixel_x_offset", 40.);
		binning.add<unsigned int>("n_bins_y", 400);
		binning.add<double>("y_min", -20.);
		binning.add<double>("y_max", 20.);

		desc.add<edm::ParameterSetDescription>("binning", binning);
	}

	descriptions.add("ppsAlignmentConfigESSource", desc);
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