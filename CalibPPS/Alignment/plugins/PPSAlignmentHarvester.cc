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
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

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
    void dqmEndJob(DQMStore::IBooker &, DQMStore::IGetter &) override {};

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
        void write(std::ostream &os) const;
        CTPPSLocalTrackLite apply(const CTPPSLocalTrackLite &tr) const;
    };

    struct AlignmentResults : public std::map<unsigned int, AlignmentResult>
    {
        void write(std::ostream &os) const;
        // int add();
        CTPPSLocalTrackLiteCollection apply(const CTPPSLocalTrackLiteCollection &input) const;
    };

    struct AlignmentResultsCollection : public std::map<std::string, AlignmentResult>
    {
        void write(std::ostream &os) const;
        // int load():
    };

    // ------------ y alignment ------------
    TF1 *ff_fit;

    double findMax();
    TGraphErrors* buildModeGraph(const TH2D *h2_y_vs_x, bool aligned, unsigned int fill, 
                                 unsigned int xangle, unsigned int rp);

    // ------------ other member data ------------
    std::string folder_;
};

// -------------------------------- Alignment classes methods --------------------------------

PPSAlignmentHarvester::AlignmentResult::AlignmentResult(double _sh_x, double _sh_x_unc, double _sh_y, 
                                                        double _sh_y_unc, double _rot_z, double _rot_z_unc)
    : sh_x(_sh_x), sh_x_unc(_sh_x_unc), sh_y(_sh_y), sh_y_unc(_sh_y_unc), rot_z(_rot_z), rot_z_unc(_rot_z_unc)
{}

void PPSAlignmentHarvester::AlignmentResult::write(std::ostream &os) const
{
    os << std::showpos << std::fixed << std::setprecision(3) << "sh_x=" << sh_x 
    << ",sh_x_unc=" << sh_x_unc << ",sh_y=" << sh_y << std::setprecision(4) << ",rot_x=" << rot_x 
    << ",rot_x_unc=" << rot_x_unc << ",rot_y=" << rot_y << ",rot_y_unc=" << rot_y_unc
    << ",rot_z=" << rot_z << ",rot_z_unc=" << rot_z_unc << "\n";
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

void PPSAlignmentHarvester::AlignmentResults::write(std::ostream &os) const
{
    for (auto &p : *this)
    {
        os << "id=" << p.first << ",";
        p.second.write(os);
    }
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

void PPSAlignmentHarvester::AlignmentResultsCollection::write(std::ostream &os) const
{
    for (auto &p : *this)
    {
        os << "\n[" << p.first.c_str() << "]\n";
        p.second.write(os);
    }
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

TGraphErrors* PPSAlignmentHarvester::buildModeGraph(const TH2D *h2_y_vs_x, bool aligned, unsigned int fill, 
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

    for (int bix = 1; bix <= h2_y_vs_x->GetNbinsX(); bix++)
    {
        const double x = h2_y_vs_x->GetXaxis()->GetBinCenter(bix);
        const double x_unc = h2_y_vs_x->GetXaxis()->GetBinWidth(bix) / 2;

        char buf[100];
        sprintf(buf, "h_y_x=%.3f", x);
		TH1D *h_y = h2_y_vs_x->ProjectionY(buf, bix, bix);

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

// -------------------------------- PPSAlignmentHarvester methods --------------------------------

PPSAlignmentHarvester::PPSAlignmentHarvester(const edm::ParameterSet &iConfig)
    : folder_(iConfig.getParameter<std::string>("folder"))
{
    // ------------ y alignment ------------
    ff_fit = new TF1("ff_fit", "[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3] + [4]*x");
}

// void PPSAlignmentHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
// {

// }

DEFINE_FWK_MODULE(PPSAlignmentHarvester);