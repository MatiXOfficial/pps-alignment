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

#include <Math/Rotation3D.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentHarvester : public DQMEDHarvester
{
public:
    PPSAlignmentHarvester(const edm::ParameterSet &);
    ~PPSAlignmentHarvester() override {};

private:
    void dqmEndJob(DQMStore::IBooker &, DQMStore::IGetter &) override;
    void dqmEndRun(DQMStore::IBooker &, DQMStore::IGetter &, edm::Run const &, edm::EventSetup const &);

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
        // int add();
        CTPPSLocalTrackLiteCollection apply(const CTPPSLocalTrackLiteCollection &input) const;
    };
    friend std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResults ar);

    struct AlignmentResultsCollection : public std::map<std::string, AlignmentResults>
    {
        // int load():
    };
    friend std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResultsCollection arc);

    // ------------ y alignment ------------
    TF1 *ff_fit;

    double findMax();
    TGraphErrors* buildModeGraph(MonitorElement *h2_y_vs_x, bool aligned, unsigned int fill, 
                                 unsigned int xangle, unsigned int rp);

    void yAlignment(DQMStore::IBooker &, DQMStore::IGetter &, const edm::EventSetup &);

    // ------------ other member data ------------
    std::string folder_;
};

// -------------------------------- Alignment classes methods --------------------------------

PPSAlignmentHarvester::AlignmentResult::AlignmentResult(double _sh_x, double _sh_x_unc, double _sh_y, 
                                                        double _sh_y_unc, double _rot_z, double _rot_z_unc)
    : sh_x(_sh_x), sh_x_unc(_sh_x_unc), sh_y(_sh_y), sh_y_unc(_sh_y_unc), rot_z(_rot_z), rot_z_unc(_rot_z_unc)
{}

std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResult ar)
{
    os << std::showpos << std::fixed << std::setprecision(3) << "sh_x=" << ar.sh_x << ",sh_x_unc=" 
    << ar.sh_x_unc << ",sh_y=" << ar.sh_y << ",sh_y_unc=" << ar.sh_y_unc << std::setprecision(4) 
    << ",rot_x=" << ar.rot_x << ",rot_x_unc=" << ar.rot_x_unc << ",rot_y=" << ar.rot_y 
    << ",rot_y_unc=" << ar.rot_y_unc << ",rot_z=" << ar.rot_z << ",rot_z_unc=" << ar.rot_z_unc << "\n";
    return os;
}

CTPPSLocalTrackLite PPSAlignmentHarvester::AlignmentResult::apply(const CTPPSLocalTrackLite &tr) const
{
    ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> v(tr.x(), tr.y(), 0.);
    ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> s(sh_x, -sh_y, 0.);
    ROOT::Math::RotationZYX R(-rot_z, rot_y, rot_x);
    v = R * v + s;

    return CTPPSLocalTrackLite(tr.rpId(), v.x(), 0., v.y(), 0., tr.tx(), tr.txUnc(), tr.ty(), tr.tyUnc(),
            tr.chiSquaredOverNDF(), tr.pixelTrackRecoInfo(), tr.numberOfPointsUsedForFit(), tr.time(), tr.timeUnc());
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

std::ostream &operator<<(std::ostream &os, PPSAlignmentHarvester::AlignmentResultsCollection arc)
{
    for (auto &p : arc)
    {
        os << "\n[" << p.first.c_str() << "]\n";
        os << p.second;
    }
    return os;
}


// -------------------------------- y alignment methods --------------------------------

double PPSAlignmentHarvester::findMax()
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

TGraphErrors* PPSAlignmentHarvester::buildModeGraph(MonitorElement *h2_y_vs_x, bool aligned, unsigned int fill, 
                                                    unsigned int xangle, unsigned int rp)
{
    std::map<unsigned int, double> mymf; // probably to be removed
    if (aligned)
    {
        mymf[23] = 3.5, mymf[3] = 4.5, mymf[103] = 5.5, mymf[123] = 4.8;
    }
    else
    {
        if (fill <= 6778)   // before TS1
		{ 
			mymf[23] = 7.2; mymf[3] = 8.3;

			if (xangle == 130) { mymf[103] = 9.0; mymf[123] = 8.0; }
			if (xangle == 160) { mymf[103] = 9.0; mymf[123] = 8.0; }
		} 
        else if (fill <= 7145)  // after TS1
		{ 
			mymf[23] = 7.5; mymf[3] = 7.8;

			if (xangle == 130) { mymf[103] = 8.2; mymf[123] = 8.0; }
			if (xangle == 160) { mymf[103] = 8.2; mymf[123] = 8.0; }
		} 
        else    // after TS2
		{ 
			mymf[23] = 7.5; mymf[3] = 7.0;

			if (xangle == 130) { mymf[103] = 7.4; mymf[123] = 8.0; }
			if (xangle == 160) { mymf[103] = 7.4; mymf[123] = 8.0; }
		}
    }

    const double yMaxFit = mymf[rp];
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

        double y_mode = findMax();
		const double y_mode_fit_unc = ff_fit->GetParameter(2) / 10;
		const double y_mode_sys_unc = 0.030;
		double y_mode_unc = std::sqrt(y_mode_fit_unc*y_mode_fit_unc + y_mode_sys_unc*y_mode_sys_unc);

		const double chiSqThreshold = (aligned) ? 1000. : 50.;

		const bool valid = ! (std::fabs(y_mode_unc) > 5. || std::fabs(y_mode) > 20. || ff_fit->GetChisquare() / ff_fit->GetNDF() > chiSqThreshold);

		if (!valid)
			continue;

		int idx = g_y_mode_vs_x->GetN();
		g_y_mode_vs_x->SetPoint(idx, x, y_mode);
		g_y_mode_vs_x->SetPointError(idx, x_unc, y_mode_unc);
    }

    return g_y_mode_vs_x;
}

void PPSAlignmentHarvester::yAlignment(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, const edm::EventSetup &iSetup)
{
    // bool useAuxFits = true;

    edm::ESHandle<PPSAlignmentConfig> cfg;
    iSetup.get<PPSAlignmentConfigRcd>().get(cfg);

    struct RPData
    {
        std::string name;
        unsigned int id;
        std::string sectorName;
        double slope;
        double sh_x;
    };

    std::vector<RPData> rpData = {
        { "L_2_F", 23,  "sector 45", (cfg->xangle() == 160) ? 0.19 : 0.17, -42. },
		{ "L_1_F",  3,  "sector 45", (cfg->xangle() == 160) ? 0.19 : 0.18, -3.6 },
		{ "R_1_F", 103, "sector 56", (cfg->xangle() == 160) ? 0.40 : 0.34, -2.8 },
		{ "R_2_F", 123, "sector 56", (cfg->xangle() == 160) ? 0.39 : 0.34, -41.9 }
    };

    // prepare results
    AlignmentResultsCollection results;

    TF1 *ff = new TF1("ff", "[0] + [1]*(x - [2])");
	TF1 *ff_sl_fix = new TF1("ff_sl_fix", "[0] + [1]*(x - [2])");

    // processing
    for (const auto &rpd : rpData)
    {
        edm::LogInfo("y_alignment") << rpd.name << "\n";
        auto *h2_y_vs_x = iGetter.get(folder_ + "/worker/" + rpd.sectorName + "/multiplicity selection/" + rpd.name + "/h2_y_vs_x");

        if (h2_y_vs_x == nullptr)
        {
            edm::LogInfo("y_alignment") << "    cannot load data, skipping\n";
            continue;
        }

        auto *g_y_cen_vs_x = buildModeGraph(h2_y_vs_x, cfg->aligned(), cfg->fill(), cfg->xangle(), rpd.id);

        if (g_y_cen_vs_x->GetN() < 5)
            continue;

        const double xMin = cfg->alignment_y_ranges()[rpd.id].x_min;
		const double xMax = cfg->alignment_y_ranges()[rpd.id].x_max;
        edm::LogInfo("y_alignment") << "    x_min = " << std::fixed << std::setprecision(3) << xMin
        << ", x_max = " << xMax << "\n";

        double sh_x = rpd.sh_x;
        double slope = rpd.slope;

        // if (useAuxFits)
        // {
        // }

        ff->SetParameters(0., 0., 0.);
        ff->FixParameter(2, -sh_x);
        ff->SetLineColor(2);
        g_y_cen_vs_x->Fit(ff, "Q", "", xMin, xMax);

        const double a = ff->GetParameter(1), a_unc = ff->GetParError(1);
		const double b = ff->GetParameter(0), b_unc = ff->GetParError(0);

        edm::LogInfo("y_alignment") << "    slope (fitted)" << std::fixed << std::setprecision(3) << a << "\n";

        results["y_alignment"][rpd.id] = AlignmentResult(0., 0., b, b_unc, 0., 0.);

        ff_sl_fix->SetParameters(0., 0., 0.);
		ff_sl_fix->FixParameter(1, slope);
		ff_sl_fix->FixParameter(2, -sh_x);
		ff_sl_fix->SetLineColor(4);
		g_y_cen_vs_x->Fit(ff_sl_fix, "Q+", "", xMin, xMax);

        const double b_fs = ff_sl_fix->GetParameter(0), b_fs_unc = ff_sl_fix->GetParError(0);

		results["y_alignment_sl_fix"][rpd.id] = AlignmentResult(0., 0., b_fs, b_fs_unc, 0., 0.);

        TGraph *g_results = new TGraph();
		g_results->SetPoint(0, sh_x, 0.);
		g_results->SetPoint(1, a, a_unc);
		g_results->SetPoint(2, b, b_unc);
		g_results->SetPoint(3, b_fs, b_fs_unc);
    }

    // write results
    edm::LogInfo("y_alignment_results") << results;
    std::cout << "y alignment results:\n" << results;
}

// -------------------------------- PPSAlignmentHarvester methods --------------------------------

PPSAlignmentHarvester::PPSAlignmentHarvester(const edm::ParameterSet &iConfig)
    : folder_(iConfig.getParameter<std::string>("folder"))
{
    // ------------ y alignment ------------
    ff_fit = new TF1("ff_fit", "[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3] + [4]*x");
}

void PPSAlignmentHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{}

void PPSAlignmentHarvester::dqmEndRun(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter, 
                                      edm::Run const &, edm::EventSetup const &iSetup)
{
    yAlignment(iBooker, iGetter, iSetup);
}


DEFINE_FWK_MODULE(PPSAlignmentHarvester);