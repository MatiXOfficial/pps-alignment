/****************************************************************************
 *
 *  CalibPPS/Alignment/plugins/PPSAlignmentHarvester.cc
 *
 *  Description : PPS Alignment DQM harvester
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

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TSystemFile.h"
#include "TSpline.h"
#include "TCanvas.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentHarvester : public DQMEDHarvester
{
public:
	PPSAlignmentHarvester(const edm::ParameterSet &iConfig);
	~PPSAlignmentHarvester() override;

private:
	void dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) override;
	void dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, edm::Run const &iRun, 
	               edm::EventSetup const &iSetup);

	// ------------ x alignment ------------
	static int fitProfile(TProfile *p, double x_mean, double x_rms, double &sl, double &sl_unc);
	static TDirectory* findDirectoryWithName(TDirectory *dir, std::string searchName);
	TGraphErrors* buildGraphFromDirectory(TDirectory *dir, const RPConfig &rpd);
	TGraphErrors* buildGraphFromMonitorElements(DQMStore::IGetter &iGetter, const RPConfig &rpd,
	                                            const std::vector<MonitorElement*> &mes);
	int doMatch(TGraphErrors *g_ref, TGraphErrors *g_test, const SelectionRange &range_ref, 
	            const SelectionRange &range_test, double sh_min, double sh_max, double sh_step,
	            double &sh_best, double &sh_best_unc);

	void xAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg,
	                const edm::ESHandle<PPSAlignmentConfig> &cfg_ref, int seqPos);

	std::map<unsigned int, double> sh_x_map;

	// ------------ x alignment relative ------------
	void xAlignmentRelative(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg, int seqPos);

	// ------------ y alignment ------------
	static double findMax(TF1 *ff_fit);
	TGraphErrors* buildModeGraph(MonitorElement *h2_y_vs_x, const edm::ESHandle<PPSAlignmentConfig> &cfg, 
	                             double x_min_mode, double x_max_mode);

	void yAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg, int seqPos);

	// ------------ other member data and methods ------------
	void debugPlots(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg);

	const std::string folder_;
	const bool debug_;
	TFile *debug_file_;
};

// -------------------------------- x alignment methods --------------------------------

int PPSAlignmentHarvester::fitProfile(TProfile *p, double x_mean, double x_rms, double &sl, double &sl_unc)
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

TDirectory* PPSAlignmentHarvester::findDirectoryWithName(TDirectory *dir, std::string searchName)
{
	TIter next(dir->GetListOfKeys());
	std::queue<TDirectory *> dirQueue;
	TObject *o;
	while ((o = next()))
	{
		TKey *k = (TKey *) o;
		
		std::string name = k->GetName();
		if (name == searchName)
			return dir;
		else if (name == "EventInfo")
			continue;

		if (((TSystemFile *) k)->IsDirectory())
			dirQueue.push((TDirectory *) k->ReadObj());
	}

	while(!dirQueue.empty())
	{
		TDirectory *resultDir = findDirectoryWithName(dirQueue.front(), searchName);
		dirQueue.pop();
		if (resultDir != nullptr)
			return resultDir;
	}

	return nullptr;
}

TGraphErrors* PPSAlignmentHarvester::buildGraphFromDirectory(TDirectory *dir, const RPConfig &rpd)
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

		TDirectory *d_slice = (TDirectory *) k->ReadObj();

		TH1D *h_y = (TH1D *) d_slice->Get("h_y");
		TProfile *p_y_diffFN_vs_y = (TProfile *) d_slice->Get("p_y_diffFN_vs_y");

		double y_cen = h_y->GetMean();
		double y_width = h_y->GetRMS();

		y_cen += rpd.y_cen_add;
		y_width *= rpd.y_width_mult;

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
	g->Sort();

	return g;
}

TGraphErrors* PPSAlignmentHarvester::buildGraphFromMonitorElements(DQMStore::IGetter &iGetter, const RPConfig &rpd,
                                                                   const std::vector<MonitorElement*> &mes)
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

			y_cen += rpd.y_cen_add;
			y_width *= rpd.y_width_mult;
			
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
	g->Sort();

	return g;
}

int PPSAlignmentHarvester::doMatch(TGraphErrors *g_ref, TGraphErrors *g_test, const SelectionRange &range_ref, 
                                   const SelectionRange &range_test, double sh_min, double sh_max, double sh_step,
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
	<< "sh_best = (" << sh_best << " +- " << sh_best_unc << ") mm";

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

void PPSAlignmentHarvester::xAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg,
                                       const edm::ESHandle<PPSAlignmentConfig> &cfg_ref, int seqPos)
{
	TDirectory *xAliDir = nullptr;
	if (debug_)
		xAliDir = debug_file_->mkdir((std::to_string(seqPos + 1) + ": x alignment").c_str());

	// prepare results
	CTPPSRPAlignmentCorrectionsData results;

	for (auto ref : cfg->matchingReferenceDatasets())
	{
		TDirectory *refDir = nullptr;
		if (debug_)
			refDir = xAliDir->mkdir(std::regex_replace(ref, std::regex("/"), "_").c_str());

		TFile *f_ref = TFile::Open(ref.c_str());
		TDirectory *ad_ref = findDirectoryWithName((TDirectory *) f_ref, cfg->sectorConfig45().name);
		if (ad_ref == nullptr)
		{
			edm::LogWarning("x_alignment") << "could not find reference dataset " << ref;
			continue;
		}
		edm::LogInfo("x_alignment") << "loading reference dataset from " << ad_ref->GetPath();

		for (const auto &sdp : { std::make_pair(cfg->sectorConfig45(), cfg_ref->sectorConfig45()), 
		                         std::make_pair(cfg->sectorConfig56(), cfg_ref->sectorConfig56()) })
		{
			const auto &sd = sdp.first;
			for (const auto &rpdp : { std::make_pair(sd.rp_F, sdp.second.rp_F), 
			                          std::make_pair(sd.rp_N, sdp.second.rp_F) })
			{
				const auto &rpd = rpdp.first;

				auto *d_ref = (TDirectory *) ad_ref->Get((sd.name + "/near_far/x slices, " + rpd.position).c_str());
				if (d_ref == nullptr)
				{
					edm::LogWarning("x_alignment") << "could not load d_ref";
					continue;
				}

				auto mes_test = iGetter.getAllContents(folder_ + "/" + sd.name + "/near_far/x slices, " + rpd.position);
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
				TGraphErrors *g_ref = buildGraphFromDirectory(d_ref, rpdp.second);

				if (debug_)
					gDirectory = rpDir->mkdir("fits_test");
				TGraphErrors *g_test = buildGraphFromMonitorElements(iGetter, rpd, mes_test);

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
				                shiftRange.x_max, cfg->x_ali_sh_step(), sh, sh_unc);
				if (r == 0)
				{
					CTPPSRPAlignmentCorrectionData rpResult(sh, sh_unc, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
					results.setRPCorrection(rpd.id, rpResult);
					edm::LogInfo("x_alignment") << std::fixed << std::setprecision(3) 
					<< "Setting sh_x of " << rpd.name << " to " << sh;
					sh_x_map[rpd.id] = sh;
				}
			}
		}
		delete f_ref;
	}

	edm::LogInfo("x_alignment_results") << seqPos + 1 << ": x_alignment:\n"<< results;
}

// -------------------------------- x alignment relative methods --------------------------------

void PPSAlignmentHarvester::xAlignmentRelative(DQMStore::IGetter &iGetter, 
                                               const edm::ESHandle<PPSAlignmentConfig> &cfg, int seqPos)
{
	TDirectory *xAliRelDir = nullptr;
	if (debug_)
		xAliRelDir = debug_file_->mkdir((std::to_string(seqPos + 1) + ": x_alignment_relative").c_str());

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
			sectorDir = xAliRelDir->mkdir(sd.name.c_str());
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

		const double sh_x_N = sh_x_map[sd.rp_N.id];
		double slope = sd.slope;

		ff->SetParameters(0., slope, 0.);
		ff->FixParameter(2, -sh_x_N);
		ff->SetLineColor(2);
		p_x_diffFN_vs_x_N->Fit(ff, "Q", "", xMin, xMax);

		const double a = ff->GetParameter(1), a_unc = ff->GetParError(1);
		const double b = ff->GetParameter(0), b_unc = ff->GetParError(0);

		edm::LogInfo("x_alignment_relative") << sd.name << ":\n" << std::fixed << std::setprecision(3)
		<< "    x_min = " << xMin << ", x_max = " << xMax << "\n"
		<< "    sh_x_N = " << sh_x_N << ", slope (fix) = " << slope << ", slope (fitted) = " << a;

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
	edm::LogInfo("x_alignment_relative_results") << seqPos + 1 << ": x_alignment_relative:\n" << results
	<< seqPos + 1 << ": x_alignment_relative_sl_fix:\n" << results_sl_fix;
}

// -------------------------------- y alignment methods --------------------------------

double PPSAlignmentHarvester::findMax(TF1 *ff_fit)
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

TGraphErrors* PPSAlignmentHarvester::buildModeGraph(MonitorElement *h2_y_vs_x, 
                                                    const edm::ESHandle<PPSAlignmentConfig> &cfg, 
                                                    double x_min_mode, double x_max_mode)
{
	TDirectory *d_top = nullptr;
	if (debug_) 
		d_top = gDirectory;

	TF1 *ff_fit = new TF1("ff_fit", "[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3] + [4]*x");

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

		double xMin = x_min_mode, xMax = x_max_mode;	// previously: xMax = yMaxFit

		h_y->Fit(ff_fit, "Q", "", xMin, xMax);

		ff_fit->ReleaseParameter(4);
		double w = std::min(4., 2. * ff_fit->GetParameter(2));
		xMin = ff_fit->GetParameter(1) - w;
		xMax = std::min(x_max_mode, ff_fit->GetParameter(1) + w);

		h_y->Fit(ff_fit, "Q", "", xMin, xMax);

		if (debug_)
			h_y->Write("h_y");

		double y_mode = findMax(ff_fit);
		const double y_mode_fit_unc = ff_fit->GetParameter(2) / 10;
		const double y_mode_sys_unc = cfg->y_mode_sys_unc();
		double y_mode_unc = std::sqrt(y_mode_fit_unc*y_mode_fit_unc + y_mode_sys_unc*y_mode_sys_unc);

		const double chiSqThreshold = cfg->chiSqThreshold();

		const bool valid = ! (std::fabs(y_mode_unc) > cfg->y_mode_unc_max_valid() 
		                      || std::fabs(y_mode) > cfg->y_mode_max_valid() 
		                      || ff_fit->GetChisquare() / ff_fit->GetNDF() > chiSqThreshold);

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

void PPSAlignmentHarvester::yAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg, 
                                       int seqPos)
{
	TDirectory *yAliDir = nullptr;
	if (debug_)
		yAliDir = debug_file_->mkdir((std::to_string(seqPos + 1) + ": y_alignment").c_str());

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
				rpDir = yAliDir->mkdir(rpd.name.c_str());
				gDirectory = rpDir->mkdir("x");
			}

			auto *h2_y_vs_x = iGetter.get(folder_ + "/" + sd.name + "/multiplicity selection/" + rpd.name + "/h2_y_vs_x");

			if (h2_y_vs_x == nullptr)
			{
				edm::LogWarning("y_alignment") << "    cannot load data for " << rpd.name << ", skipping";
				continue;
			}

			auto *g_y_cen_vs_x = buildModeGraph(h2_y_vs_x, cfg, rpd.x_min_mode, rpd.x_max_mode);

			if (g_y_cen_vs_x->GetN() < 5)
			{
				edm::LogInfo("y_alignment") << "    insufficient data, skipping";
				continue;
			}

			const double xMin = cfg->alignment_y_ranges()[rpd.id].x_min;
			const double xMax = cfg->alignment_y_ranges()[rpd.id].x_max;

			const double sh_x = sh_x_map[rpd.id];
			double slope = rpd.slope;

			ff->SetParameters(0., 0., 0.);
			ff->FixParameter(2, -sh_x);
			ff->SetLineColor(2);
			g_y_cen_vs_x->Fit(ff, "Q", "", xMin, xMax);

			const double a = ff->GetParameter(1), a_unc = ff->GetParError(1);
			const double b = ff->GetParameter(0), b_unc = ff->GetParError(0);

			edm::LogInfo("y_alignment") << rpd.name  << ":\n" << std::fixed << std::setprecision(3)
			<< "    x_min = " << xMin << ", x_max = " << xMax << "\n"
			<< "    sh_x = " << sh_x << ", slope (fix) = " << slope << ", slope (fitted) = " << a;

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
	edm::LogInfo("y_alignment_results") << seqPos + 1 << ": y_alignment:\n" << results 
	<< seqPos + 1 << ": y_alignment_sl_fix:\n" << results_sl_fix;
}

// -------------------------------- PPSAlignmentHarvester methods --------------------------------

void PPSAlignmentHarvester::debugPlots(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg)
{
	TDirectory *debugPlotsDir = debug_file_->mkdir("misc");
	for (const auto &sd : { cfg->sectorConfig45(), cfg->sectorConfig56() })
	{
		for (const auto &rpd : { sd.rp_F, sd.rp_N })
		{
			gDirectory = debugPlotsDir->mkdir(rpd.name.c_str());

			auto *tmp = iGetter.get(folder_ + "/" + sd.name + "/profiles/" + rpd.name + "/m_p_y_vs_x_aft_sel");
			if (tmp == nullptr)
				edm::LogWarning("debug_plots") << "could not load the m_p_y_vs_x_aft_sel Tprofile";
			auto *profile = tmp->getTProfile();

			int bins = profile->GetXaxis()->GetNbins();
			double xMin = profile->GetXaxis()->GetXmin();
			double xMax = profile->GetXaxis()->GetXmax();

			auto *h_entries = new TH1D("h_entries", ";x", bins, xMin, xMax);
			auto *h_mean = new TH1D("h_mean", ";x", bins, xMin, xMax);
			auto *h_stddev = new TH1D("h_stddev", ";x", bins, xMin, xMax);
			
			for (int bin = 1; bin < bins; bin++)
			{
				int N = profile->GetBinEntries(bin);
				h_entries->SetBinContent(bin, N);

				h_mean->SetBinContent(bin, profile->GetBinContent(bin));
				h_mean->SetBinError(bin, profile->GetBinError(bin));

				profile->SetErrorOption("s");
				double s = profile->GetBinError(bin);
				h_stddev->SetBinContent(bin, s);
				h_stddev->SetBinError(bin, (N > 0) ? s / std::sqrt(2. * N) : 0.);

				profile->SetErrorOption();
			}

			h_entries->Write("h_entries");
			h_mean->Write("h_mean");
			h_stddev->Write("h_stddev");
		}
	}
}

PPSAlignmentHarvester::PPSAlignmentHarvester(const edm::ParameterSet &iConfig)
	: folder_(iConfig.getParameter<std::string>("folder")),
	  debug_(iConfig.getParameter<bool>("debug"))
{}

PPSAlignmentHarvester::~PPSAlignmentHarvester()
{}

void PPSAlignmentHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{}

void PPSAlignmentHarvester::dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                      edm::Run const &iRun, edm::EventSetup const &iSetup)
{
	edm::ESHandle<PPSAlignmentConfig> cfg;
	iSetup.get<PPSAlignmentConfigRcd>().get(cfg);
	
	edm::ESHandle<PPSAlignmentConfig> cfg_ref;
	iSetup.get<PPSAlignmentConfigRcd>().get("reference", cfg_ref);
	
	if (debug_)
	{
		debug_file_ = new TFile(("debug_" + std::to_string(iRun.run()) + ".root").c_str(), "recreate");
		debugPlots(iGetter, cfg);
	}

	// setting default sh_x values from config
	for (const auto sd : { cfg->sectorConfig45(), cfg->sectorConfig56() })
	{
		for (const auto rpd : { sd.rp_N, sd.rp_F })
		{
			edm::LogInfo("harvester_dqmEndRun") << std::fixed << std::setprecision(3) 
			<< "Setting sh_x of " << rpd.name << " to " << rpd.sh_x;
			sh_x_map[rpd.id] = rpd.sh_x;
		}
	}

	for (unsigned int i = 0; i < cfg->sequence().size(); i++)
	{
		if (cfg->sequence()[i] == "x alignment")
			xAlignment(iGetter, cfg, cfg_ref, i);
		else if (cfg->sequence()[i] == "x alignment relative")
			xAlignmentRelative(iGetter, cfg, i);
		else if (cfg->sequence()[i] == "y alignment")
			yAlignment(iGetter, cfg, i);
		else
			edm::LogError("harvester_dqmEndRun") << cfg->sequence()[i] << " is a wrong method name.";
	}

	if (debug_)
		delete debug_file_;
}

DEFINE_FWK_MODULE(PPSAlignmentHarvester);