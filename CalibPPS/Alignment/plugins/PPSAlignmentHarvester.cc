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

// -------------------------------- PPSAlignmentHarvester methods --------------------------------

PPSAlignmentHarvester::PPSAlignmentHarvester(const edm::ParameterSet &iConfig) 
{
    // ------------ y alignment ------------
    ff_fit = new TF1("ff_fit", "[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3] + [4]*x");
}

void PPSAlignmentHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter)
{
    // ------------ y alignment ------------
}

DEFINE_FWK_MODULE(PPSAlignmentHarvester);