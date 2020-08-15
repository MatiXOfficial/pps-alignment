/****************************************************************************
 *
 *  CalibPPS/Alignment/plugins/PPSAlignmentXHarvester.cc
 *
 *  Description : PPS Alignment DQM harvester - x alignment
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

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <regex>

#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TSpline.h"
#include "TCanvas.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentXHarvester : public DQMEDHarvester
{
public:
	PPSAlignmentXHarvester(const edm::ParameterSet &iConfig);
	~PPSAlignmentXHarvester() override;

private:
	void dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) override;
	void dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, edm::Run const &, 
	               edm::EventSetup const &iSetup) override;

	// ------------ x alignment ------------
	static int fitProfile(TProfile *p, double x_mean, double x_rms, double &sl, double &sl_unc);
	TGraphErrors* buildGraphFromDirectory(TDirectory *dir, bool aligned, unsigned int rpId);
	TGraphErrors* buildGraphFromMonitorElements(DQMStore::IGetter &iGetter, 
	                                            const std::vector<MonitorElement*> &mes, 
	                                            bool aligned, unsigned int rpId);
	int doMatch(TGraphErrors *g_ref, TGraphErrors *g_test, const SelectionRange &range_ref, 
	            const SelectionRange &range_test, double sh_min, double sh_max, double &sh_best, 
	            double &sh_best_unc);

	void xAlignment(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
	                const edm::ESHandle<PPSAlignmentConfig> &cfg, 
	                const edm::ESHandle<PPSAlignmentConfig> &cfg_ref);

	// ------------ other member data and methods ------------
	const std::string folder_;
	const bool debug_;
};

// -------------------------------- x alignment methods --------------------------------

int PPSAlignmentXHarvester::fitProfile(TProfile *p, double x_mean, double x_rms, double &sl, double &sl_unc)
{
	if (p->GetEntries() < 50)
		return 1;

	unsigned int n_reasonable = 0;
	for (int bi = 1; bi <= p->GetNbinsX(); bi++)
	{
		if (p->GetBinEntries(bi) < 5)
		{
			p->SetBinContent(bi, 0.);
			p->SetBinError(bi, 0.);
		} 
		else 
		{
			n_reasonable++;
		}
	}

	if (n_reasonable < 10)
		return 2;

	double xMin = x_mean - x_rms, xMax = x_mean + x_rms;

	TF1 *ff_pol1 = new TF1("ff_pol1", "[0] + [1]*x");

	ff_pol1->SetParameter(0., 0.);
	p->Fit(ff_pol1, "Q", "", xMin, xMax);

	sl = ff_pol1->GetParameter(1);
	sl_unc = ff_pol1->GetParError(1);

	return 0;
}

TGraphErrors* PPSAlignmentXHarvester::buildGraphFromDirectory(TDirectory *dir, bool aligned, unsigned int rpId)
{
	TGraphErrors *g = new TGraphErrors();

	TIter next(dir->GetListOfKeys());
	TObject *o;
	while ((o = next()))
	{
		TKey *k = (TKey *) o;

		std::string name = k->GetName();
		size_t d = name.find("-");
		const double x_min = std::stod(name.substr(0, d));
		const double x_max = std::stod(name.substr(d+1));

		//printf("  %s, %.3f, %.3f\n", name.c_str(), x_min, x_max);

		TDirectory *d_slice = (TDirectory *) k->ReadObj();

		TH1D *h_y = (TH1D *) d_slice->Get("h_y");
		TProfile *p_y_diffFN_vs_y = (TProfile *) d_slice->Get("p_y_diffFN_vs_y");

		double y_cen = h_y->GetMean();
		double y_width = h_y->GetRMS();

		if (aligned)
		{
			y_cen += ((rpId < 100) ? -0.2 : -0.4);
		}
		else 
		{
			y_cen += ((rpId < 100) ? -0.3 : -0.8);
			y_width *= ((rpId < 100) ? 1.1 : 1.0);
		}

		double sl=0., sl_unc=0.;
		int fr = fitProfile(p_y_diffFN_vs_y, y_cen, y_width, sl, sl_unc);
		if (fr != 0)
			continue;

		if (debug_)
			p_y_diffFN_vs_y->Write(name.c_str());

		int idx = g->GetN();
		g->SetPoint(idx, (x_max + x_min)/2., sl);
		g->SetPointError(idx, (x_max - x_min)/2., sl_unc);
	}

	return g;
}

TGraphErrors* PPSAlignmentXHarvester::buildGraphFromMonitorElements(DQMStore::IGetter &iGetter, 
                                                                   const std::vector<MonitorElement*> &mes, 
                                                                   bool aligned, unsigned int rpId)
{
	TGraphErrors *g = new TGraphErrors();

	for (auto *me : mes)
	{
		if (me->getName() == "h_y")
		{
			std::string parentPath = me->getPathname();
			size_t parentPos = parentPath.substr(0, parentPath.size() - 1).find_last_of("/") + 1;
			std::string parentName = parentPath.substr(parentPos);
			size_t d = parentName.find("-");
			const double xMin = std::stod(parentName.substr(0, d));
			const double xMax = std::stod(parentName.substr(d + 1));

			TH1D *h_y = me->getTH1D();

			auto *p_y_diffFN_vs_y_monitor = iGetter.get(parentPath + "p_y_diffFN_vs_y");
			if (p_y_diffFN_vs_y_monitor == nullptr)
			{
				edm::LogWarning("x_alignment") << "could not find p_y_diffFN_vs_y in: " << parentPath;
				continue;
			}
			TProfile *p_y_diffFN_vs_y = p_y_diffFN_vs_y_monitor->getTProfile();

			double y_cen = h_y->GetMean();
			double y_width = h_y->GetRMS();

			if (aligned)
			{
				y_cen += (rpId < 100) ? -0.2 : -0.4;
			}   
			else
			{
				y_cen += (rpId < 100) ? -0.3 : -0.8;
				y_width *= (rpId < 100) ? 1.1 : 1.0;
			}
			
			double sl = 0., sl_unc = 0.;
			int fr = fitProfile(p_y_diffFN_vs_y, y_cen, y_width, sl, sl_unc);
			if (fr != 0)
				continue;

			if (debug_)
				p_y_diffFN_vs_y->Write(parentName.c_str());

			int idx = g->GetN();
			g->SetPoint(idx, (xMax + xMin) / 2., sl);
			g->SetPointError(idx, (xMax - xMin) / 2., sl_unc);
		}
	}

	return g;
}

int PPSAlignmentXHarvester::doMatch(TGraphErrors *g_ref, TGraphErrors *g_test, const SelectionRange &range_ref, 
                                   const SelectionRange &range_test, double sh_min, double sh_max, 
                                   double &sh_best, double &sh_best_unc)
{
	// require minimal number of points
	if (g_ref->GetN() < 5 || g_test->GetN() < 5)
		return 1;

	// print config
	edm::LogInfo("x_alignment") << std::fixed << std::setprecision(3) 
	<< "ref: x_min = " << range_ref.x_min << ", x_max = " << range_ref.x_max << "\n"
	<< "test: x_min = " << range_test.x_min << ", x_max = " << range_test.x_max;

	// make spline from g_ref
	TSpline3 *s_ref = new TSpline3("s_ref", g_ref->GetX(), g_ref->GetY(), g_ref->GetN());

	// book match-quality graphs
	TGraph *g_n_points = new TGraph(); g_n_points->SetName("g_n_points"); g_n_points->SetTitle(";sh;N");
	TGraph *g_chi_sq = new TGraph(); g_chi_sq->SetName("g_chi_sq"); g_chi_sq->SetTitle(";sh;S2");
	TGraph *g_chi_sq_norm = new TGraph(); g_chi_sq_norm->SetName("g_chi_sq_norm"); g_chi_sq_norm->SetTitle(";sh;S2 / N");

	// optimalisation variables
	double S2_norm_best = 1E100;

	double sh_step = 0.010;	// mm
	for (double sh = sh_min; sh <= sh_max; sh += sh_step)
	{
		// calculate chi^2
		int n_points = 0;
		double S2 = 0.;

		for (int i = 0; i < g_test->GetN(); ++i)
		{
			const double x_test = g_test->GetX()[i];
			const double y_test = g_test->GetY()[i];
			const double y_test_unc = g_test->GetErrorY(i);

			const double x_ref = x_test + sh;

			if (x_ref < range_ref.x_min || x_ref > range_ref.x_max || x_test < range_test.x_min || x_test > range_test.x_max)
				continue;

			const double y_ref = s_ref->Eval(x_ref);

			int js = -1, jg = -1;
			double xs = -1E100, xg = +1E100;
			for (int j = 0; j < g_ref->GetN(); ++j)
			{
				const double x = g_ref->GetX()[j];
				if (x < x_ref && x > xs)
				{
					xs = x;
					js = j;
				}
				if (x > x_ref && x < xg)
				{
					xg = x;
					jg = j;
				}
			}
			if (jg == -1)
				jg = js;

			const double y_ref_unc = ( g_ref->GetErrorY(js) + g_ref->GetErrorY(jg) ) / 2.;

			n_points++;
			const double S2_inc = pow(y_test - y_ref, 2.) / (y_ref_unc*y_ref_unc + y_test_unc*y_test_unc);
			S2 += S2_inc;
		}

		// update best result
		double S2_norm = S2 / n_points;

		if (S2_norm < S2_norm_best)
		{
			S2_norm_best = S2_norm;
			sh_best = sh;
		}

		// fill in graphs
		int idx = g_n_points->GetN();
		g_n_points->SetPoint(idx, sh, n_points);
		g_chi_sq->SetPoint(idx, sh, S2);
		g_chi_sq_norm->SetPoint(idx, sh, S2_norm);
	}

	TF1 *ff_pol2 = new TF1("ff_pol2", "[0] + [1]*x + [2]*x*x");

	// determine uncertainty
	double fit_range = 0.5;	// mm
	g_chi_sq->Fit(ff_pol2, "Q", "", sh_best - fit_range, sh_best + fit_range);
	sh_best_unc = 1. / sqrt(ff_pol2->GetParameter(2));

	// print results
	edm::LogInfo("x_alignment") << std::fixed << std::setprecision(3) 
	<< "sh_best = (" << sh_best << " +- " << sh_best_unc << " mm";

	if (debug_)
	{
		// save graphs
		g_n_points->Write();
		g_chi_sq->Write();
		g_chi_sq_norm->Write();

		// save results
		TGraph *g_results = new TGraph();
		g_results->SetName("g_results");
		g_results->SetPoint(0, sh_best, sh_best_unc);
		g_results->SetPoint(1, range_ref.x_min, range_ref.x_max);
		g_results->SetPoint(2, range_test.x_min, range_test.x_max);
		g_results->Write();

		// save debug canvas
		TGraphErrors *g_test_shifted = new TGraphErrors(*g_test);
		for (int i = 0; i < g_test_shifted->GetN(); ++i)
		{
			g_test_shifted->GetX()[i] += sh_best;
		}

		TCanvas *c_cmp = new TCanvas("c_cmp");
		g_ref->SetLineColor(1);
		g_ref->SetName("g_ref");
		g_ref->Draw("apl");

		g_test->SetLineColor(4);
		g_test->SetName("g_test");
		g_test->Draw("pl");

		g_test_shifted->SetLineColor(2);
		g_test_shifted->SetName("g_test_shifted");
		g_test_shifted->Draw("pl");
		c_cmp->Write();

		delete c_cmp;
	}

	// clean up
	delete s_ref;

	return 0;
}

void PPSAlignmentXHarvester::xAlignment(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                        const edm::ESHandle<PPSAlignmentConfig> &cfg,
                                        const edm::ESHandle<PPSAlignmentConfig> &cfg_ref)
{
	TFile *debugFile = nullptr;
	if (debug_)
		debugFile = new TFile("x_alignment_debug.root", "recreate");

	// prepare results
	CTPPSRPAlignmentCorrectionsData results;

	for (auto ref : cfg->matchingReferenceDatasets())
	{
		edm::LogInfo("x_alignment") << "reference dataset: " << ref;

		TDirectory *refDir = nullptr;
		if (debug_)
			refDir = debugFile->mkdir(std::regex_replace(ref, std::regex("/"), "_").c_str());

		TFile *f_ref = TFile::Open(ref.c_str());

		for (const auto &sd : { cfg->sectorConfig45(), cfg->sectorConfig56() })
		{
			for (const auto &rpd : { sd.rp_F, sd.rp_N })
			{
				auto *d_ref = (TDirectory *) f_ref->Get((sd.name + "/near_far/x slices, " + rpd.position).c_str());
				auto mes_test = iGetter.getAllContents(folder_ + "/" + sd.name + "/near_far/x slices, " + rpd.position);

				if (d_ref == nullptr)
				{
					edm::LogWarning("x_alignment") << "could not load d_ref";
					continue;
				}

				if (mes_test.empty())
				{
					edm::LogWarning("x_alignment") << "could not load mes_test";
					continue;
				}

				TDirectory *rpDir = nullptr;
				if (debug_)
				{
					rpDir = refDir->mkdir(rpd.name.c_str());
					gDirectory = rpDir->mkdir("fits_ref");
				}
				TGraphErrors *g_ref = buildGraphFromDirectory(d_ref, cfg_ref->aligned(), rpd.id);

				if (debug_)
					gDirectory = rpDir->mkdir("fits_test");
				TGraphErrors *g_test = buildGraphFromMonitorElements(iGetter, mes_test, cfg->aligned(), rpd.id);

				if (debug_)
				{
					gDirectory = rpDir;
					g_ref->Write("g_ref");
					g_test->Write("g_test");
				}

				const auto &shiftRange = cfg->matchingShiftRanges()[rpd.id];
				double sh = 0., sh_unc = 0.;
				int r = doMatch(g_ref, g_test, cfg_ref->alignment_x_meth_o_ranges()[rpd.id], 
								cfg->alignment_x_meth_o_ranges()[rpd.id], shiftRange.x_min, 
								shiftRange.x_max, sh, sh_unc);
				if (r == 0)
				{
					CTPPSRPAlignmentCorrectionData rpResult(sh, sh_unc, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
					results.setRPCorrection(rpd.id, rpResult);
				}
				
				iBooker.setCurrentFolder(folder_ + "/" + sd.name + "/sh_x/" + rpd.name);
				auto *sh_me = iBooker.bookFloat("sh_x");
				sh_me->Fill(sh);
			}
		}
		delete f_ref;
	}

	edm::LogInfo("x_alignment_results") << "x_alignment_meth_o:\n"<< results;
	
	if (debug_)
		delete debugFile;
}



// -------------------------------- PPSAlignmentXHarvester methods --------------------------------

PPSAlignmentXHarvester::PPSAlignmentXHarvester(const edm::ParameterSet &iConfig)
	: folder_(iConfig.getParameter<std::string>("folder")),
	  debug_(iConfig.getParameter<bool>("debug"))
{
	std::cout << "!! Rozpoczeto x !!\n";
}

PPSAlignmentXHarvester::~PPSAlignmentXHarvester()
{}

void PPSAlignmentXHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{}

void PPSAlignmentXHarvester::dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                      edm::Run const &, edm::EventSetup const &iSetup)
{
	edm::ESHandle<PPSAlignmentConfig> cfg;
	iSetup.get<PPSAlignmentConfigRcd>().get(cfg);
	
	edm::ESHandle<PPSAlignmentConfig> cfg_ref;
	iSetup.get<PPSAlignmentConfigRcd>().get("reference", cfg_ref);

	xAlignment(iBooker, iGetter, cfg, cfg_ref);
	std::cout << "!! Zakonczono x !!\n";
}

DEFINE_FWK_MODULE(PPSAlignmentXHarvester);