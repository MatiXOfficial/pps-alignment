/****************************************************************************
 *
 *  CalibPPS/Alignment/plugins/PPSAlignmentYHarvester.cc
 *
 *  Description : PPS Alignment DQM harvester - y alignment
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
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentYHarvester : public DQMEDHarvester
{
public:
	PPSAlignmentYHarvester(const edm::ParameterSet &iConfig);
	~PPSAlignmentYHarvester() override;

private:
	void dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) override;
	void dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, edm::Run const &, 
	               edm::EventSetup const &iSetup) override;

	// ------------ y alignment ------------
	static double findMax(TF1 *ff_fit);
	TGraphErrors* buildModeGraph(MonitorElement *h2_y_vs_x, bool aligned, double _yMaxFit);

	void yAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg);

	// ------------ other member data and methods ------------
	const std::string folder_;
	const bool debug_;
};

// -------------------------------- y alignment methods --------------------------------

double PPSAlignmentYHarvester::findMax(TF1 *ff_fit)
{
	const double mu = ff_fit->GetParameter(1);
	const double si = ff_fit->GetParameter(2);

	// unreasonable fit?
	if (si > 25. || std::fabs(mu) > 100.)
		return 1E100;

	double xMax = 1E100;
	double yMax = -1E100;
	for (double x = mu - si; x <= mu + si; x += 0.001)
	{
		double y = ff_fit->Eval(x);
		if (y > yMax)
		{
			xMax = x;
			yMax = y;
		}
	}

	return xMax;
}

TGraphErrors* PPSAlignmentYHarvester::buildModeGraph(MonitorElement *h2_y_vs_x, bool aligned, double _yMaxFit)
{
	TDirectory *d_top = nullptr;
	if (debug_) 
		d_top = gDirectory;

	TF1 *ff_fit = new TF1("ff_fit", "[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3] + [4]*x");

	const double yMaxFit = _yMaxFit;
	TGraphErrors *g_y_mode_vs_x = new TGraphErrors();

	for (int bix = 1; bix <= h2_y_vs_x->getNbinsX(); bix++)
	{
		const double x = h2_y_vs_x->getTH2D()->GetXaxis()->GetBinCenter(bix);
		const double x_unc = h2_y_vs_x->getTH2D()->GetXaxis()->GetBinWidth(bix) / 2;

		char buf[100];
		sprintf(buf, "h_y_x=%.3f", x);
		TH1D *h_y = h2_y_vs_x->getTH2D()->ProjectionY(buf, bix, bix);

		if (h_y->GetEntries() < 300)
			continue;

		if (debug_)
		{
			sprintf(buf, "x=%.3f", x);
			gDirectory = d_top->mkdir(buf);
		}

		double conMax = -1.;
		double conMax_x = 0.;
		for (int biy = 1; biy < h_y->GetNbinsX(); biy++)
		{
			if (h_y->GetBinContent(biy) > conMax)
			{
				conMax = h_y->GetBinContent(biy);
				conMax_x = h_y->GetBinCenter(biy);
			}
		}

		ff_fit->SetParameters(conMax, conMax_x, h_y->GetRMS() * 0.75, 0., 0.);
		ff_fit->FixParameter(4, 0.);

		double xMin = 2., xMax = yMaxFit;
		if (aligned)
			xMin = -2., xMax = +3.;

		h_y->Fit(ff_fit, "Q", "", xMin, xMax);

		ff_fit->ReleaseParameter(4);
		double w = std::min(4., 2. * ff_fit->GetParameter(2));
		xMin = ff_fit->GetParameter(1) - w;
		xMax = std::min(yMaxFit, ff_fit->GetParameter(1) + w);

		h_y->Fit(ff_fit, "Q", "", xMin, xMax);

		if (debug_)
			h_y->Write("h_y");

		double y_mode = findMax(ff_fit);
		const double y_mode_fit_unc = ff_fit->GetParameter(2) / 10;
		const double y_mode_sys_unc = 0.030;
		double y_mode_unc = std::sqrt(y_mode_fit_unc*y_mode_fit_unc + y_mode_sys_unc*y_mode_sys_unc);

		const double chiSqThreshold = (aligned) ? 1000. : 50.;

		const bool valid = ! (std::fabs(y_mode_unc) > 5. || std::fabs(y_mode) > 20. || ff_fit->GetChisquare() / ff_fit->GetNDF() > chiSqThreshold);

		if (debug_)
		{
			TGraph *g_data = new TGraph();
			g_data->SetPoint(0, y_mode, y_mode_unc);
			g_data->SetPoint(1, ff_fit->GetChisquare(), ff_fit->GetNDF());
			g_data->SetPoint(2, valid, 0.);
			g_data->Write("g_data");
		}

		if (!valid)
			continue;

		int idx = g_y_mode_vs_x->GetN();
		g_y_mode_vs_x->SetPoint(idx, x, y_mode);
		g_y_mode_vs_x->SetPointError(idx, x_unc, y_mode_unc);
	}

	return g_y_mode_vs_x;
}

void PPSAlignmentYHarvester::yAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg)
{
	TFile *debugFile = nullptr;
	if (debug_)
		debugFile = new TFile("y_alignment_debug.root", "recreate");

	// prepare results
	CTPPSRPAlignmentCorrectionsData results;
	CTPPSRPAlignmentCorrectionsData results_sl_fix;

	TF1 *ff = new TF1("ff", "[0] + [1]*(x - [2])");
	TF1 *ff_sl_fix = new TF1("ff_sl_fix", "[0] + [1]*(x - [2])");

	// processing
	for (const auto &sd : { cfg->sectorConfig45(), cfg->sectorConfig56() })
	{
		for (const auto &rpd : { sd.rp_F, sd.rp_N })
		{
			TDirectory *rpDir = nullptr;
			if (debug_)
			{
				rpDir = debugFile->mkdir(rpd.name.c_str());
				gDirectory = rpDir->mkdir("x");
			}
			
			auto *sh_x_me = iGetter.get(folder_ + "/" + sd.name + "/sh_x/" + rpd.name + "/sh_x");
			if (sh_x_me == nullptr)
				std::cout << "NIE ZNALEZIONO sh_x dla " << rpd.name << "\n============================\n";
			else
				std::cout << "JEST!!!!! dla " << rpd.name << " i = " << sh_x_me->getFloatValue() << "\n++++++++++++++++++++++++++++++++\n";

			auto *h2_y_vs_x = iGetter.get(folder_ + "/" + sd.name + "/multiplicity selection/" + rpd.name + "/h2_y_vs_x");

			if (h2_y_vs_x == nullptr)
			{
				edm::LogWarning("y_alignment") << "    cannot load data for " << rpd.name << ", skipping";
				continue;
			}

			auto *g_y_cen_vs_x = buildModeGraph(h2_y_vs_x, cfg->aligned(), cfg->yMaxFit()[rpd.id]);

			if (g_y_cen_vs_x->GetN() < 5)
			{
				edm::LogInfo("y_alignment") << "    insufficient data, skipping";
				continue;
			}

			const double xMin = cfg->alignment_y_ranges()[rpd.id].x_min;
			const double xMax = cfg->alignment_y_ranges()[rpd.id].x_max;

			double sh_x = rpd.sh_x;
			double slope = rpd.slope;

			ff->SetParameters(0., 0., 0.);
			ff->FixParameter(2, -sh_x);
			ff->SetLineColor(2);
			g_y_cen_vs_x->Fit(ff, "Q", "", xMin, xMax);

			const double a = ff->GetParameter(1), a_unc = ff->GetParError(1);
			const double b = ff->GetParameter(0), b_unc = ff->GetParError(0);

			edm::LogInfo("y_alignment") << rpd.name << std::fixed << std::setprecision(3) << ":\n"
			<< "    x_min = " << xMin << ", x_max = " << xMax << "\n"
			<< "    sh_x = " << sh_x << ", slope (fix) = " << slope << "\n"
			<< "    slope (fitted) = " << a;

			CTPPSRPAlignmentCorrectionData rpResult(0., 0., b, b_unc, 0., 0., 0., 0., 0., 0., 0., 0.);
			results.setRPCorrection(rpd.id, rpResult);

			ff_sl_fix->SetParameters(0., 0., 0.);
			ff_sl_fix->FixParameter(1, slope);
			ff_sl_fix->FixParameter(2, -sh_x);
			ff_sl_fix->SetLineColor(4);
			g_y_cen_vs_x->Fit(ff_sl_fix, "Q+", "", xMin, xMax);

			const double b_fs = ff_sl_fix->GetParameter(0), b_fs_unc = ff_sl_fix->GetParError(0);
			
			CTPPSRPAlignmentCorrectionData rpResult_sl_fix(0., 0., b_fs, b_fs_unc, 0., 0., 0., 0., 0., 0., 0., 0.);
			results_sl_fix.setRPCorrection(rpd.id, rpResult_sl_fix);

			if (debug_)
			{
				gDirectory = rpDir;

				g_y_cen_vs_x->Write("g_y_cen_vs_x");

				TGraph *g_results = new TGraph();
				g_results->SetPoint(0, sh_x, 0.);
				g_results->SetPoint(1, a, a_unc);
				g_results->SetPoint(2, b, b_unc);
				g_results->SetPoint(3, b_fs, b_fs_unc);
				g_results->Write("g_results");
			}
		}
	}

	// write results
	edm::LogInfo("y_alignment_results") << "y_alignment:\n" << results << "y_alignment_sl_fix:\n" << results_sl_fix;
	
	if (debug_)
		delete debugFile;
}

// -------------------------------- PPSAlignmentYHarvester methods --------------------------------

PPSAlignmentYHarvester::PPSAlignmentYHarvester(const edm::ParameterSet &iConfig)
	: folder_(iConfig.getParameter<std::string>("folder")),
	  debug_(iConfig.getParameter<bool>("debug"))
{
	std::cout << "!! Rozpoczeto y !!\n";
}

PPSAlignmentYHarvester::~PPSAlignmentYHarvester()
{}

void PPSAlignmentYHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{}

void PPSAlignmentYHarvester::dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                      edm::Run const &, edm::EventSetup const &iSetup)
{
	edm::ESHandle<PPSAlignmentConfig> cfg;
	iSetup.get<PPSAlignmentConfigRcd>().get(cfg);

	yAlignment(iGetter, cfg);
	std::cout << "!! Zakonczono y !!\n";
}

DEFINE_FWK_MODULE(PPSAlignmentYHarvester);