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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <regex>

#include <Math/Rotation3D.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
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
	void dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, edm::Run const &, 
	               edm::EventSetup const &iSetup);

	// ------------ structures ------------
	struct AlignmentResult
	{
		double sh_x = 0., sh_x_unc = 0.;    // mm
		double sh_y = 0., sh_y_unc = 0.;    // mm

		double rot_x = 0., rot_x_unc = 0.;  // rad
		double rot_y = 0., rot_y_unc = 0.;  // rad
		double rot_z = 0., rot_z_unc = 0.;  // rad

		AlignmentResult(double _sh_x = 0., double _sh_x_unc = 0., double _sh_y = 0., 
		                double _sh_y_unc = 0., double _rot_z = 0., double _rot_z_unc = 0.);

		CTPPSLocalTrackLite apply(const CTPPSLocalTrackLite &tr) const;
	};
	friend std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResult ar);

	struct AlignmentResults : public std::map<unsigned int, AlignmentResult>
	{
		CTPPSLocalTrackLiteCollection apply(const CTPPSLocalTrackLiteCollection &input) const;
	};
	friend std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResults ar);

	struct AlignmentResultsCollection : public std::map<std::string, AlignmentResults>
	{};
	friend std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResultsCollection arc);

	// ------------ x alignment ------------
	static int fitProfile(TProfile *p, double x_mean, double x_rms, double &sl, double &sl_unc);
	TGraphErrors* buildGraphFromDirectory(TDirectory *dir, bool aligned, unsigned int rpId);
	TGraphErrors* buildGraphFromMonitorElements(DQMStore::IGetter &iGetter, 
	                                            const std::vector<MonitorElement*> &mes, 
	                                            bool aligned, unsigned int rpId);
	int doMatch(TGraphErrors *g_ref, TGraphErrors *g_test, const SelectionRange &range_ref, 
	            const SelectionRange &range_test, double sh_min, double sh_max, double &sh_best, 
	            double &sh_best_unc);

	void xAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg,
	                const edm::ESHandle<PPSAlignmentConfig> &cfg_ref);

	// ------------ x alignment relative ------------
	void xAlignmentRelative(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg);

	// ------------ y alignment ------------
	static double findMax(TF1 *ff_fit);
	TGraphErrors* buildModeGraph(MonitorElement *h2_y_vs_x, bool aligned, double _yMaxFit);

	void yAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg);

	// ------------ other member data and methods ------------
	void debugPlots(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
	                const edm::ESHandle<PPSAlignmentConfig> &cfg);

	const std::string folder_;
	const bool debug_;
	TFile *debug_file_;
};

// -------------------------------- Alignment classes methods --------------------------------

PPSAlignmentHarvester::AlignmentResult::AlignmentResult(double _sh_x, double _sh_x_unc, double _sh_y, 
                                                        double _sh_y_unc, double _rot_z, double _rot_z_unc)
	: sh_x(_sh_x), sh_x_unc(_sh_x_unc), sh_y(_sh_y), sh_y_unc(_sh_y_unc), rot_z(_rot_z), rot_z_unc(_rot_z_unc)
{}

CTPPSLocalTrackLite PPSAlignmentHarvester::AlignmentResult::apply(const CTPPSLocalTrackLite &tr) const
{
	ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> v(tr.x(), tr.y(), 0.);
	ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> s(sh_x, -sh_y, 0.);
	ROOT::Math::RotationZYX R(-rot_z, rot_y, rot_x);
	v = R * v + s;

	return CTPPSLocalTrackLite(tr.rpId(), v.x(), 0., v.y(), 0., tr.tx(), tr.txUnc(), tr.ty(), tr.tyUnc(),
	                           tr.chiSquaredOverNDF(), tr.pixelTrackRecoInfo(), tr.numberOfPointsUsedForFit(), 
	                           tr.time(), tr.timeUnc());
}

std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResult ar)
{
	os << std::showpos << std::fixed << std::setprecision(3) << "sh_x=" << ar.sh_x << ",sh_x_unc=" 
	<< ar.sh_x_unc << ",sh_y=" << ar.sh_y << ",sh_y_unc=" << ar.sh_y_unc << std::setprecision(4) 
	<< ",rot_x=" << ar.rot_x << ",rot_x_unc=" << ar.rot_x_unc << ",rot_y=" << ar.rot_y 
	<< ",rot_y_unc=" << ar.rot_y_unc << ",rot_z=" << ar.rot_z << ",rot_z_unc=" << ar.rot_z_unc << "\n";
	return os;
}

CTPPSLocalTrackLiteCollection PPSAlignmentHarvester::AlignmentResults::apply(const CTPPSLocalTrackLiteCollection &input) const
{
	CTPPSLocalTrackLiteCollection output;
	for (auto &t : input)
	{
		auto ait = find(t.rpId());
		if (ait == end())
			throw cms::Exception("alignment") << "No alignment data for RP " << t.rpId();
		
		output.emplace_back(ait->second.apply(t));
	}
	return output;
}

std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResults ar)
{
	for (auto &p : ar)
	{
		os << "id=" << p.first << ",";
		os << p.second;
	}
	return os;
}

std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResultsCollection arc)
{
	for (auto &p : arc)
	{
		os << "\n[" << p.first.c_str() << "]\n";
		os << p.second;
	}
	return os;
}

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

TGraphErrors* PPSAlignmentHarvester::buildGraphFromDirectory(TDirectory *dir, bool aligned, unsigned int rpId)
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

TGraphErrors* PPSAlignmentHarvester::buildGraphFromMonitorElements(DQMStore::IGetter &iGetter, 
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

int PPSAlignmentHarvester::doMatch(TGraphErrors *g_ref, TGraphErrors *g_test, const SelectionRange &range_ref, 
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

void PPSAlignmentHarvester::xAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg,
                                       const edm::ESHandle<PPSAlignmentConfig> &cfg_ref)
{
	TDirectory *xAliDir = nullptr;
	if (debug_)
		xAliDir = debug_file_->mkdir("x alignment");

	// prepare results
	AlignmentResultsCollection results;

	for (auto ref : cfg->matchingReferenceDatasets())
	{
		edm::LogInfo("x_alignment") << "reference dataset: " << ref;

		TDirectory *refDir = nullptr;
		if (debug_)
			refDir = xAliDir->mkdir(std::regex_replace(ref, std::regex("/"), "_").c_str());

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
					results["x_alignment_meth_o"][rpd.id] = AlignmentResult(sh, sh_unc);
			}
		}
		delete f_ref;
	}

	edm::LogInfo("x_alignment_results") << results;
}

// -------------------------------- x alignment relative methods --------------------------------

void PPSAlignmentHarvester::xAlignmentRelative(DQMStore::IGetter &iGetter, 
                                               const edm::ESHandle<PPSAlignmentConfig> &cfg)
{
	TDirectory *xAliRelDir = nullptr;
	if (debug_)
		xAliRelDir = debug_file_->mkdir("x_alignment_relative");

	// prepare results
	AlignmentResultsCollection results;

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
			edm::LogInfo("x_alignment_relative") << "    insufficient data, skipping";
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

		results["x_alignment_relative"][sd.rp_N.id] = AlignmentResult(+b/2., b_unc/2., 0., 0., 0., 0.);
		results["x_alignment_relative"][sd.rp_F.id] = AlignmentResult(-b/2., b_unc/2., 0., 0., 0., 0.);

		ff_sl_fix->SetParameters(0., slope, 0.);
		ff_sl_fix->FixParameter(1, slope);
		ff_sl_fix->FixParameter(2, -sh_x_N);
		ff_sl_fix->SetLineColor(4);
		p_x_diffFN_vs_x_N->Fit(ff_sl_fix, "Q+", "", xMin, xMax);

		const double b_fs = ff_sl_fix->GetParameter(0), b_fs_unc = ff_sl_fix->GetParError(0);

		results["x_alignment_relative_sl_fix"][sd.rp_N.id] = AlignmentResult(+b_fs/2., b_fs_unc/2., 0., 0., 0., 0.);
		results["x_alignment_relative_sl_fix"][sd.rp_F.id] = AlignmentResult(-b_fs/2., b_fs_unc/2., 0., 0., 0., 0.);

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
	edm::LogInfo("x_alignment_relative_results") << results;
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

TGraphErrors* PPSAlignmentHarvester::buildModeGraph(MonitorElement *h2_y_vs_x, bool aligned, double _yMaxFit)
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

void PPSAlignmentHarvester::yAlignment(DQMStore::IGetter &iGetter, const edm::ESHandle<PPSAlignmentConfig> &cfg)
{
	TDirectory *yAliDir = nullptr;
	if (debug_)
		yAliDir = debug_file_->mkdir("y_alignment");

	// prepare results
	AlignmentResultsCollection results;

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

			results["y_alignment"][rpd.id] = AlignmentResult(0., 0., b, b_unc, 0., 0.);

			ff_sl_fix->SetParameters(0., 0., 0.);
			ff_sl_fix->FixParameter(1, slope);
			ff_sl_fix->FixParameter(2, -sh_x);
			ff_sl_fix->SetLineColor(4);
			g_y_cen_vs_x->Fit(ff_sl_fix, "Q+", "", xMin, xMax);

			const double b_fs = ff_sl_fix->GetParameter(0), b_fs_unc = ff_sl_fix->GetParError(0);

			results["y_alignment_sl_fix"][rpd.id] = AlignmentResult(0., 0., b_fs, b_fs_unc, 0., 0.);

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
	edm::LogInfo("y_alignment_results") << results;
}

// -------------------------------- PPSAlignmentHarvester methods --------------------------------

void PPSAlignmentHarvester::debugPlots(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                       const edm::ESHandle<PPSAlignmentConfig> &cfg)
{
	for (const auto &sd : { cfg->sectorConfig45(), cfg->sectorConfig56() })
	{
		for (const auto &rpd : { sd.rp_F, sd.rp_N })
		{
			auto *tmp = iGetter.get(folder_ + "/" + sd.name + "/profiles/" + rpd.name + "/m_p_y_vs_x_aft_sel");
			if (tmp == nullptr)
				edm::LogWarning("debug_plots") << "could not load the m_p_y_vs_x_aft_sel Tprofile";
			auto *profile = tmp->getTProfile();

			int bins = profile->GetXaxis()->GetNbins();
			double xMin = profile->GetXaxis()->GetXmin();
			double xMax = profile->GetXaxis()->GetXmax();

			iBooker.setCurrentFolder(folder_ + "/" + sd.name + "/profiles/" + rpd.name);
			auto *h_entries = iBooker.book1DD("h_entries", ";x", bins, xMin, xMax);
			auto *h_mean = iBooker.book1DD("h_mean", ";x", bins, xMin, xMax);
			auto *h_stddev = iBooker.book1DD("h_stddev", ";x", bins, xMin, xMax);
			
			for (int bin = 1; bin < bins; bin++)
			{
				int N = profile->GetBinEntries(bin);
				h_entries->setBinContent(bin, N);

				h_mean->setBinContent(bin, profile->GetBinContent(bin));
				h_mean->setBinError(bin, profile->GetBinError(bin));

				profile->SetErrorOption("s");
				double s = profile->GetBinError(bin);
				h_stddev->setBinContent(bin, s);
				h_stddev->setBinError(bin, (N > 0) ? s / std::sqrt(2. * N) : 0.);

				profile->SetErrorOption();
			}
		}
	}
}

PPSAlignmentHarvester::PPSAlignmentHarvester(const edm::ParameterSet &iConfig)
	: folder_(iConfig.getParameter<std::string>("folder")),
	  debug_(iConfig.getParameter<bool>("debug"))
{
	if (debug_)
		debug_file_ = new TFile("debug.root", "recreate");
}

PPSAlignmentHarvester::~PPSAlignmentHarvester()
{
	if (debug_)
		delete debug_file_;
}

void PPSAlignmentHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{}

void PPSAlignmentHarvester::dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                      edm::Run const &, edm::EventSetup const &iSetup)
{
	edm::ESHandle<PPSAlignmentConfig> cfg;
	iSetup.get<PPSAlignmentConfigRcd>().get(cfg);
	
	edm::ESHandle<PPSAlignmentConfig> cfg_ref;
	iSetup.get<PPSAlignmentConfigRcd>().get("reference", cfg_ref);
	
	if (debug_)
		debugPlots(iBooker, iGetter, cfg);

	xAlignment(iGetter, cfg, cfg_ref);
	xAlignmentRelative(iGetter, cfg);
	yAlignment(iGetter, cfg);
}

DEFINE_FWK_MODULE(PPSAlignmentHarvester);