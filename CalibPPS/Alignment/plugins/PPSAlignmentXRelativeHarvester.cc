/****************************************************************************
 *
 *  CalibPPS/Alignment/plugins/PPSAlignmentXRelativeHarvester.cc
 *
 *  Description : PPS Alignment DQM harvester - x relative alignment
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/PPSObjects/interface/CTPPSRPAlignmentCorrectionData.h"
#include "CondFormats/PPSObjects/interface/CTPPSRPAlignmentCorrectionsData.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <string>
#include <iostream>
#include <iomanip>

#include "TGraphErrors.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentXRelativeHarvester : public DQMEDHarvester
{
public:
	PPSAlignmentXRelativeHarvester(const edm::ParameterSet &iConfig);
	~PPSAlignmentXRelativeHarvester() override;

private:
	void dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) override;
	void dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, edm::Run const &, 
	               edm::EventSetup const &iSetup) override;

	// ------------ x alignment relative ------------
	void xAlignmentRelative(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg);

	// ------------ other member data and methods ------------
	void debugPlots(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
	                const edm::ESHandle<PPSAlignmentConfig> &cfg);

	const std::string folder_;
	const bool debug_;
};

// -------------------------------- x alignment relative methods --------------------------------

void PPSAlignmentXRelativeHarvester::xAlignmentRelative(DQMStore::IGetter &iGetter, 
                                               const edm::ESHandle<PPSAlignmentConfig> &cfg)
{
	TFile *debugFile = nullptr;
	if (debug_)
		debugFile = new TFile("x_alignment_relative_debug.root", "recreate");

	// prepare results
	CTPPSRPAlignmentCorrectionsData results;
	CTPPSRPAlignmentCorrectionsData results_sl_fix;

	TF1 *ff = new TF1("ff", "[0] + [1]*(x - [2])");
	TF1 *ff_sl_fix = new TF1("ff_sl_fix", "[0] + [1]*(x - [2])");

	// processing
	for (const auto &sd : { cfg->sectorConfig45(), cfg->sectorConfig56() })
	{
		TDirectory *sectorDir = nullptr;
		if (debug_)
		{
			sectorDir = debugFile->mkdir(sd.name.c_str());
			gDirectory = sectorDir;
		}
		
		auto *p_x_diffFN_vs_x_N_monitor = iGetter.get(folder_ + "/" + sd.name + "/near_far/p_x_diffFN_vs_x_N");
		if (p_x_diffFN_vs_x_N_monitor == nullptr)
		{
			edm::LogWarning("x_alignment_relative") << "    cannot load data, skipping";
			continue;
		}
		TProfile *p_x_diffFN_vs_x_N = p_x_diffFN_vs_x_N_monitor->getTProfile();

		if (p_x_diffFN_vs_x_N->GetEntries() < 100)
		{
			edm::LogInfo("x_alignment_relative") << "insufficient data, skipping";
			continue;
		}

		const double xMin = cfg->alignment_x_relative_ranges()[sd.rp_N.id].x_min;
		const double xMax = cfg->alignment_x_relative_ranges()[sd.rp_N.id].x_max;

		edm::LogInfo("x_alignment_relative") << sd.name << std::fixed << std::setprecision(3) << ":\n"
		<< "    x_min = " << xMin << ", x_max = " << xMax;

		double slope = sd.slope;
		double sh_x_N = sd.rp_N.sh_x;

		ff->SetParameters(0., slope, 0.);
		ff->FixParameter(2, -sh_x_N);
		ff->SetLineColor(2);
		p_x_diffFN_vs_x_N->Fit(ff, "Q", "", xMin, xMax);

		const double a = ff->GetParameter(1), a_unc = ff->GetParError(1);
		const double b = ff->GetParameter(0), b_unc = ff->GetParError(0);

		CTPPSRPAlignmentCorrectionData rpResult_N(+b/2., b_unc/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
		results.setRPCorrection(sd.rp_N.id, rpResult_N);
		CTPPSRPAlignmentCorrectionData rpResult_F(-b/2., b_unc/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
		results.setRPCorrection(sd.rp_F.id, rpResult_F);

		ff_sl_fix->SetParameters(0., slope, 0.);
		ff_sl_fix->FixParameter(1, slope);
		ff_sl_fix->FixParameter(2, -sh_x_N);
		ff_sl_fix->SetLineColor(4);
		p_x_diffFN_vs_x_N->Fit(ff_sl_fix, "Q+", "", xMin, xMax);

		const double b_fs = ff_sl_fix->GetParameter(0), b_fs_unc = ff_sl_fix->GetParError(0);
		
		CTPPSRPAlignmentCorrectionData rpResult_sl_fix_N(+b_fs/2., b_fs_unc/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
		results_sl_fix.setRPCorrection(sd.rp_N.id, rpResult_sl_fix_N);
		CTPPSRPAlignmentCorrectionData rpResult_sl_fix_F(-b_fs/2., b_fs_unc/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
		results_sl_fix.setRPCorrection(sd.rp_F.id, rpResult_sl_fix_F);

		if (debug_)
		{
			p_x_diffFN_vs_x_N->Write("p_x_diffFN_vs_x_N");

			TGraph *g_results = new TGraph();
			g_results->SetPoint(0, sh_x_N, 0.);
			g_results->SetPoint(1, a, a_unc);
			g_results->SetPoint(2, b, b_unc);
			g_results->SetPoint(3, b_fs, b_fs_unc);
			g_results->Write("g_results");
		}
	}

	// write results
	edm::LogInfo("x_alignment_relative_results") << "x_alignment_relative:\n" << results
	<< "x_alignment_relative_sl_fix:\n" << results_sl_fix;
	
	if (debug_)
		delete debugFile;
}

// -------------------------------- PPSAlignmentXRelativeHarvester methods --------------------------------

PPSAlignmentXRelativeHarvester::PPSAlignmentXRelativeHarvester(const edm::ParameterSet &iConfig)
	: folder_(iConfig.getParameter<std::string>("folder")),
	  debug_(iConfig.getParameter<bool>("debug"))
{
	std::cout << "!! Rozpoczeto x relative !!\n";
}

PPSAlignmentXRelativeHarvester::~PPSAlignmentXRelativeHarvester()
{}

void PPSAlignmentXRelativeHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{}

void PPSAlignmentXRelativeHarvester::dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                      edm::Run const &, edm::EventSetup const &iSetup)
{
	edm::ESHandle<PPSAlignmentConfig> cfg;
	iSetup.get<PPSAlignmentConfigRcd>().get(cfg);

	xAlignmentRelative(iGetter, cfg);
	std::cout << "!! Zakonczono x relative !!\n";
}

DEFINE_FWK_MODULE(PPSAlignmentXRelativeHarvester);